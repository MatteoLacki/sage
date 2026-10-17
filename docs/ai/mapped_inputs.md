# Peptide table and positional predictions mapped in place (2026-10-17)

Plan: necromerge2 `plans/mapped_peptide_table.md`. Everything SAGE reads per
peptide row is written once, by the step that produces it, in the layout SAGE uses,
as mmappet datasets; SAGE maps them read-only, so they are clean file pages the OS
can drop and re-read instead of anonymous memory. Row numbers (`PeptideIx`,
`peptide_row`) do not change.

## Stored compact table (`--compact-peptides`, sliced search only)

`dump_peptides --compact <dir>` builds the run's `CompactPeptides` table
(`CompactPeptides::build`: direct for non-specific digests, per length for
enzymatic ones) and stores it with `sage_cloudpath::compact_store::write`. With
`--compact` alone it skips sage's own digest and writes no parquet; with both, it
checks that the table's rows and the digest's peptides have the same
`dumped_peptides_sha256`.

`<dir>` is an mmappet dataset over the rows, exactly `CompactRows`' columns
(`compact/rows.rs`) in native (little-endian) byte order:

| `schema.txt` | column |
|---|---|
| `uint16 mass_residual` | mass bit pattern minus its 1,024-row chunk's base |
| `uint32 start` | proteome offset of the target sequence |
| `uint8 length` | |
| `uint8 flags` | decoy, semi-enzymatic, position, missed cleavages |
| `uint32 var_sites` | the two `VarSite`s' four bytes (`#[repr(C)]`), opaque outside SAGE |

plus sidecars `mass_chunks.bin` (per chunk: base bit pattern, start in the wide
masses or `u32::MAX`, two `u32`), `wide_masses.bin` (full bit patterns of chunks
spanning more than 2^16) and `meta.json` (format and version, row count, rows per
chunk, `dumped_peptides_sha256`, the digest settings, the distinct variable-mod
masses as bit patterns, a SHA-256 of the proteome).

The proteome is not stored. `compact_store::open` builds it from the run's own FASTA
(tens of MB, needed anyway by `materialize`/`resolve_proteins`) and refuses the table
unless three things match: the digest settings (decoy_tag, enzyme, generate_decoys,
max_variable_mods, peptide min/max mass, static and variable mods: the keys of the
pipeline's `dump_peptides_config.json`, compared after printing both through JSON so
`f32`s compare as the same text), the variable-mod masses, and the proteome hash
(boundaries, accessions, residues). The FASTA path itself does not matter.

In memory, a `CompactRows` column is `Owned(Vec<T>)` (built) or `Stored(SharedBytes)`
(any `AsRef<[u8]>` kept alive by the column, here a `memmap2::Mmap`), checked once for
length and alignment; every accessor goes through `as_slice`, so the search code is
the same either way. sage-core gains no dependency (mapping and JSON live in
sage-cloudpath). The stored fingerprint (`CompactPeptides::stored_fingerprint`)
replaces hashing every row's mass when positional predictions are checked.

The sliced search's memory budget is half of `MemAvailable`, which counts reclaimable
page cache: with the table mapped, slices get larger (F9477 class I 3 → 2 slices,
class II 15 → 13), so part of the saved anonymous memory is spent on bigger slices
rather than showing up as a lower peak.

### Measured, F9477 immunopeptidomics (16 threads, grid layout)

Warm page cache, nothing else running, same binary; necromerge2
`temp/mapped_table/timing*`.

| | class I built | class I mapped | class II built | class II mapped |
|---|---:|---:|---:|---:|
| table | 7.7 s build | 0.1 s map | 33.0 s build | 0.1 s map |
| anonymous at "database ready" | 1.83 GB | 0.04 GB | 7.49 GB | 0.04 GB |
| slices | 3 | 2 | 15 | 13 |
| planning | 0.6 s | 3.8 s | 2.2 s | 8.5 s |
| sliced search | 85.5 s | 87.4 s | 324.9 s | 323.6 s |
| anonymous after the search | 11.57 GB | 9.95 GB | 17.01 GB | 10.64 GB |
| peak anonymous (polled) | 21.72 GB | 20.54 GB | 25.38 GB | 20.45 GB |
| wall | 107.6 s | 103.0 s | 386.7 s | 357.2 s |

Storing the table once: `dump_peptides --compact` 10.2 s and 3.6 GB peak for class I,
43.8 s and 17.0 GB for class II (the direct build's own peak).

Planning is slower mapped (+3 s, +6 s): it reads every row's mass and length through
page faults, partly from disk when the table was evicted (0.8 and 3.3 GB read). An
`madvise(MADV_SEQUENTIAL)` on the columns is the obvious next step if it matters.

Results identical to the pages baselines (PSM key, matched fragments, pin) for both
classes, mapped and built, and on the small checks (pmsms sliced, G23 tryptic).

## Positional predictions as mmappet (`--predicted-rt`, `--predicted-iim`, `--predicted-fragment-intensity-index`)

The three files `git/featureprediction` writes are now mmappet datasets with a
`meta.json` holding `dumped_peptides_sha256` (it used to be parquet file metadata):
`rt` (float64), `iim` (float64, dense over row x charge), `start`/`end` (int64, `-1`
for no cache entry). Same values and order as the parquets; float64 stays because
`git/featureprediction`'s precursor correction fits on them.

`predicted_properties::MappedColumn<T>` maps one column (dtype checked against
`schema.txt`, clones share the mapping). `Runner` keeps the mapped columns;
`Scorer`/`RtSlabbing` get `sage_core::predicted::PredictedValues` (an `&[f64]` read
as `f32`, the cast the parquet reader made once at load) and `PredictedRanges` (the
two pointer columns; a negative start is no entry). The `Option` around RT and IIM
values is gone: the loader wrapped every value in `Some`, so a value is always
present, NaN included, exactly as before (`RtSlabbing` still sends non-finite values
to the no-prediction slab, and eviction still drops them). Nothing counts the
fragment pointers with an entry any more at load: that touched every page.

At ~620M peptides and 3 charges this is ~5 GB (RT, `Option<f32>` per row), ~15 GB
(IIM, `Option<f32>` per row and charge) and ~45 GB (pointers, 24 B `Option<(u64,
u64)>` per row and charge) of anonymous memory that become mapped pages (2.5, 15 and
30 GB of files).

### Measured, F9477 best final search (tryptic, 6.3M peptides, 18.9M pointer slots)

The featureprediction steps rerun into mmappet reproduce the existing nodes bit for bit
(RT values and fingerprint, RT tolerance, fragment start/end, RT-corrected precursors).
The new binary mapping them and the previous binary reading the parquets give
identical outputs (`temp/mapped_table/stage2.sh`, `timing.sh`):

| | parquet, loaded | mmappet, mapped |
|---|---:|---:|
| anonymous at "database and predictions ready" | 4.74 GB | 4.20 GB |
| anonymous after the search | 5.56 GB | 4.99 GB |
| peak anonymous (polled, after the predictions are released) | 6.84 GB | 6.84 GB |
| search | 121.3 s | 116.2 s |
| wall | 152 s | 144 s |

At this size the saving is ~0.55 GB; it scales with peptides x charges.
