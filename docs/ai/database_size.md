# `database_size`: how big the index gets, per peptide length (2026-10-13)

`crates/sage-cli/src/bin/database_size.rs` takes the same `--fasta` / `--config`
(`database` section JSON) / `--unimod-db-path` as `dump_peptides`, and runs sage's
own `Parameters::digest` + `build_from_peptides` once per peptide length
(`enzyme.min_len = enzyme.max_len = L`). It prints a TSV with one row per length and
a `total` row:
`length, target_peptides, decoy_peptides, fragments, pages, fragment_index_gb,
peptides_gb, compact_gb, digest_s, build_s, peak_rss_gb`. `compact_gb` is the same
peptides as 16 B `CompactPeptide` rows (`compact_peptides.md`).

## Why per length is exact

Dedup in `reorder_peptides` only merges peptides with identical sequences. A decoy is
dropped only when its sequence equals a target's, and reversal preserves length. So
peptides of different lengths never interact, and the per-length rows sum to the
full database. Verified: the F9477 tryptic config
(`nodes/write_dump_peptides_config/*/dump_peptides_config.json`,
`~/fastas/human.fasta`) totals 6,300,882 peptides and 198,932,246 fragments. That is
exactly what `run_sage` logs for the same database.

Peak memory is therefore that of the largest single length, which is what lets this
measure databases that don't fit when built whole. `peak_rss_gb` is `VmHWM`, reset
before each length by writing `5` to `/proc/self/clear_refs` (Linux only). It
includes the FASTA, which is held throughout.

Byte columns are estimates:
- **Fragments:** f32 m/z + u32 peptide index + one f32 per page.
- **Peptides:** `size_of::<Peptide>()` + sequence `Arc` (refcounts + bytes) +
  modification `Vec<f32>` capacity + protein-pointer `Vec` capacity. Protein names
  are shared with the FASTA and not counted.
- **Not counted:** allocator rounding and per-allocation overhead, so real RSS is
  higher.

## Unspecific human, immunopeptidomics mods (2026-10-13)

Config: `cleave_at = ""`, no static mods, variable `M` UNIMOD:35, `[` UNIMOD:1,
`C` UNIMOD:312, `max_variable_mods = 1`, mass 600–4000, decoys. FASTA:
`~/fastas/human.fasta` (11,460,518 residues).

| range | peptides | fragments | index + `Peptide` | digest+build | max per-length peak |
|---|---:|---:|---:|---:|---:|
| 8–12 (class I) | 149M | 2.11G | 16.9 + 25.5 = 42.4 GB | 131 s | 21.1 GB |
| 7–12 | 174M | 2.31G | 48.0 GB | 152 s | 21.1 GB |
| 8–14 | 216M | 3.51G | 66.2 GB | 197 s | 24.3 GB |
| 7–24 | 629M | 16.7G | 261 GB | 677 s | 48.4 GB |

Length 25 aborted under `ulimit -v 50000000`. Extrapolating, class II (9–25) is
about 275 GB.

What the numbers say:
- **A `Peptide` costs 161–181 B (8–12 aa)**, more than its fragments: 80–144 B at
  8 B each, 10–18 fragments with `min_ion_index = 2`.
- **A single length's build peaks at about 2x its finished size** (length 12:
  10.4 GB finished, 21.1 GB peak).
- **Saving the index would not help memory.** The build is cheap next to its size
  (tryptic: 11.6 s in `run_sage`). See `plans/sliced_search.md` in the necromerge2
  root for the layout this motivated.

`digest_s` runs at about 400% CPU on 16 cores. Two parts of sage's digest are
single-threaded:
- `group_digests` (`enzyme.rs`): a serial sort of `String`-keyed digests, then a
  grouping loop.
- `reorder_peptides`' `dedup_by`.

The tool also frees each length's database serially before the next one starts.
That alternation between serial and parallel phases is why CPU usage jumps around
in htop.
