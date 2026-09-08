# `--predicted-rt`/`--predicted-iim`/`--predicted-fragment-intensity-index`: positional, not `sequence`-keyed (2026-09-08)

`Runner::new` used to load these three files into `String`-keyed
`HashMap`s, then `run()`'s `resolve_predicted_rt`/`resolve_predicted_iim`/
`resolve_predicted_fragment_intensity` called `Peptide::to_string()` once
per peptide in the digested database (6.3M times for a real F9477 run) to
build lookup keys — `to_string()` does a `unimod::lookup_reverse` per
modified residue, real measured overhead (see the `+37%/152s` note that
used to live on this code before the RT/IIM split). Loading the files
themselves was also slow: `read_fragment_intensity_index`'s row-oriented
parquet read heap-allocated a `String` per row, 16.87M times.

## The fix: peptide row position is already a valid key

`dump_peptides.rs` and a real search's own `database.build(fasta)` call
the *literal same* `Parameters::digest` (`database.rs:162`), which already
runs `reorder_peptides` (mass-sort + dedup) — see `docs/ai/dump_peptides.md`.
Row `i` of a `dumped_peptides`/`peptides.parquet` is therefore guaranteed
to be assigned `PeptideIx(i)` by any SAGE run whose `--fasta` and
`database` config trace back to that same `dumped_peptides` node — and
necroflow's own node-identity hashing already enforces exactly that
equality (`necromerge2`'s `pipelines.py`, `dump_peptides_config` is a
direct slice of the same `cfg.sage.database` dict rendered into the real
search's own `sage_config.json`).

So `git/featureprediction`'s exporters can write `predicted_rt.parquet`/
`predicted_iim.parquet`/`fragment_intensity_for_sage.parquet` with **no**
`sequence`/`charge` key columns at all — just the values, in a fixed,
predictable row order:

- `predicted_rt.parquet`: single column `rt` (f64), row `i` = peptide row
  `i`.
- `predicted_iim.parquet`: single column `iim` (f64), dense over
  `(peptide_row, charge)` in exactly `Scorer::iim_dense_slot`'s own
  `peptide_idx * charge_span + (charge - min_charge)` layout — this
  already was `predict_iim`'s natural `np.repeat`/`np.tile` cross-join
  order, so no reshaping was needed on the Python side, only removing the
  now-redundant key columns.
- `fragment_intensity_for_sage.parquet`: `start`/`end` (i64, sentinel
  `-1` where a `(peptide_row, charge)` pair has no cache entry — this one
  genuinely can be sparse, unlike RT/IIM which always predict every
  peptide), same dense slot layout.

Reading became a straight columnar/positional load (`predicted_properties.rs`,
`fragment_intensity_cache.rs`) instead of a `HashMap` build, and
`resolve_predicted_rt`/`resolve_predicted_iim`/
`resolve_predicted_fragment_intensity` (`runner.rs`) collapsed to a length
check plus a direct copy — no per-peptide string formatting or lookup at
all.

## Safety net: a fingerprint, not a silent assumption

Trusting row position instead of a self-describing key removes the
"lookup miss" failure mode a config-drifted file would hit before (loud,
if inconvenient). To keep that a loud failure rather than silent
misattribution, every positional file carries a `dumped_peptides_sha256`
parquet file-metadata entry: a SHA-256 over the source `dumped_peptides`'
`monoisotopic` column, row order, little-endian `f32` bytes
(`dumped_peptides_fingerprint.rs`/`.py` — the two implementations were
verified to produce bit-identical hashes for the same input). `run()`
recomputes the same fingerprint over `database.peptides` and hard-errors
on any mismatch, before trusting a single row position. Deliberately
cheap: no string formatting on the SAGE side, `monoisotopic` is a plain
field on `Peptide` and a plain column in `dumped_peptides`.

Row count is checked too (`values.len()` against `db.peptides.len()` or
`db.peptides.len() * charge_span`), catching a mismatch before the
fingerprint pass even needs to run.

## Known limitation: chunked prefiltering

`prefilter_peptides`'s per-fasta-chunk digest is a strict subset of any
`dumped_peptides` node's peptide list, in no relation to its row order —
positional files can't be resolved against a chunk-scoped digest at all.
`predicted_rt`/`predicted_iim` are now unconditionally skipped during
prefiltering (same existing precedent as fragment-intensity, which this
code path already skipped regardless of the real search's config) with a
one-time warning if configured. Not exercised by this repo's actual job
configs today (`database.prefilter` isn't set `true` anywhere in
`necromerge2`'s `jobs/*.toml`) — revisit if that changes.

## Verified, real F9477 data (2026-09-08)

Regenerated `predicted_rt.parquet` (6,300,882 rows) and
`fragment_intensity_for_sage.parquet` (18,902,646 dense rows, 16,874,172
populated / 2,028,474 sentinel) with the new code against the same
`dumped_peptides`/cache inputs as an existing production run, then:

- Reindexed the old string-keyed files by `peptide_row` and diffed against
  the new positional files directly — **exact match**, 0.0 max diff on
  `rt`, all populated `(start, end)` pairs identical, every non-populated
  slot correctly `-1`.
- Reran the exact real production `run_sage` command (full F9477, 716,614
  precursors) with the new binary against the new files: **363,156 PSMs**,
  identical `(scannr, peptide, label)` id set to the stored baseline, all
  54 numeric result columns byte-identical including `hyperscore` and
  `ms2_entropy_similarity`.

Cross-language fingerprint algorithm equivalence (Python `hashlib.sha256`
vs Rust `sha2::Sha256`, same little-endian `f32` bytes) verified directly
with a small fixture before the full-scale check.
