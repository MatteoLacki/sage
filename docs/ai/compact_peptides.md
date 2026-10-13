# `CompactPeptides`: 16 B peptide rows in sage's own order (2026-10-13)

`crates/sage/src/compact.rs` is Stage 1 of `plans/sliced_search.md` (necromerge2 root):
a peptide database small enough to hold whole for unspecific immunopeptidomics
digests. Nothing in the search uses it yet.

## Layout

- **`Proteome`:** every target protein's residues concatenated (`u8`), plus
  `protein_starts: Vec<u32>` and the accessions. Human: 11.5 MB.
- **`CompactPeptide`:** 16 B (`f32` mass, `u32` proteome offset of the target
  sequence, `u8` length, `u8` flags, two `(u8 position, u8 mass id)` variable-mod
  sites, padded from 14 to 16 by alignment).
  - The flags hold decoy, semi-enzymatic, `Position` and missed cleavages (4 bits).
  - Decoys are not stored as sequences. A decoy row points at its target window and
    is read through `Peptide::reverse`, which keeps both ends and flips the interior.
  - Static mods are not stored; they are re-applied from the sequence and `Position`.
- **Row order:** exactly `Parameters::digest`'s (`reorder_peptides`: mass `total_cmp`,
  then `Peptide::initial_sort`). Row `i` is the `PeptideIx(i)` a monolithic search
  assigns, so positional files (`dumped_peptides_positional_predictions.md`) line up.

`materialize(row)` rebuilds the `Peptide` through sage's own code path, so it matches
bit for bit, mass included:
1. `Digest` → `Peptide::try_from`.
2. `apply_site` for the stored variable sites.
3. `static_mods`.
4. `monoisotopic += modification_mass()`.
5. `reverse` for decoys.

`proteins` is left empty. `resolve_proteins(&mut [Peptide])` fills it, meant for the
few peptides that end up in PSMs. It runs sage's own `EnzymeParameters::digest` per
protein and keeps the proteins whose digest yields the peptide's target form with the
same modifications. That reproduces both `group_digests`' grouping and
`reorder_peptides`' merge across positions.

## Building

`CompactPeptides::build` takes four steps:
1. **Digest each length** with sage's own `Parameters::digest`
   (`Parameters::for_single_length`). This is exact for the reason in
   `database_size.md`: lengths never interact.
2. **Compress each peptide** by diffing it against its static-mods-only
   materialization. The slots that differ are its variable-mod sites.
3. **Merge the per-length runs.** Each run is already in sage's order. Between runs
   the lengths differ, and sequences of different lengths are never equal, so mass
   then sequence bytes settle every comparison. Modifications are never needed.
   Merges run in parallel over mass ranges. The pivots are sampled from every run,
   because runs of different lengths cover different masses. Pivots taken from the
   longest run (length 25, all heavy) put nearly every row in the first chunk: 74 s
   for 673M rows, against 1.96 s for 174M rows with sampled pivots.
4. **Locate each target sequence.** A `DashMap` of distinct sequences is filled in
   parallel, then frozen with `into_read_only()`, and the proteome windows are
   scanned in parallel for the first occurrence (smallest offset, so the rows are
   deterministic). The freeze matters: `DashMap::get` takes a shard lock, and with
   11.5M window lookups per length it tripled tryptic compress time (2.7 s → 6.4 s).

Speed: profiled on length 8 of the immunopeptidomics config (27.6M peptides).
- **compress:** 10.7 s → 4.1 s, from two changes. The single-threaded `FnvHashMap`
  fill took 6.05 s; the parallel `DashMap` takes 1.23 s. Dropping the decoys' target
  sequences took 1.60 s on one thread; dropping them in parallel takes 0.31 s.
  Lengths 7–12 now compress in 4.2–4.9 s each, down from 10.6–20.1 s.
- **Freeing the digest's `Vec<Peptide>`:** 14.1 s single-threaded, 4.9 s with
  `drop_in_parallel`. glibc serialises frees of memory allocated by other threads,
  which also caps the parallel drop at about 4x.
- **What's left is sage's own `Parameters::digest`:** about 20 s per length at 3x
  parallel. The single-threaded parts are `group_digests`' sort and grouping and
  `reorder_peptides`' `dedup_by`.

`group_digests` (`enzyme.rs`) now returns early on an empty input. It used to panic
on `digests[0]`, which a single-length digest can hit (no peptide of that length).

Restrictions, all hard errors from `CompactPeptides::new` or `compress`:
- `generate_decoys: true`
- at most 2 variable mods per peptide when any are configured
- `max_len` ≤ 254
- ≤ 15 missed cleavages
- ≤ 254 distinct variable-mod masses
- residues and rows fit `u32`

## Duplicate attributes are arbitrary in sage itself

When `reorder_peptides` merges equal peptides from different digests (same sequence,
mods and mass, from different proteins or protein positions), the survivor's
`position`, `semi_enzymatic` and `missed_cleavages` belong to whichever duplicate
`par_sort_unstable_by` left first. That depends on the whole array, so a per-length
digest and a full digest can keep different ones.

- **`position`:** read only while applying mods, which is already done, so it
  doesn't matter.
- **`missed_cleavages`:** an LDA feature. It is a function of the sequence for
  fully-enzymatic and nonspecific digests. Only semi-enzymatic digests can merge
  duplicates that disagree.

The unit tests compare everything else exactly, and require each compact row's
attributes to belong to one of its duplicates.

## Verification

- **Unit tests:** `cargo test -p sage-core --lib compact`, comparing against
  `Parameters::digest` row for row on a small FASTA. Configs:
  - tryptic
  - semi-tryptic
  - nonspecific
  - nonspecific with static peptide-N-term, variable protein-N-term, peptide-C-term,
    `^Q` and two `M` masses
  - immunopeptidomics (cysteinylation, max 1 variable mod)
- **`verify_compact_peptides` binary:** per length, every row materializes back to
  sage's `Peptide` and every N-th row resolves sage's protein list. The merged order
  is checked between adjacent rows. With `--full`, the whole table is compared to a
  whole-database digest.
- **F9477 tryptic** (`dump_peptides_config.json`, `~/fastas/human.fasta`, `--full`,
  proteins checked on every row): 6,300,882 rows, 0 mismatches, sage's exact order,
  `position` differing on 195 merged peptides, `missed_cleavages` and
  `semi_enzymatic` on none. 0.10 GB against 1.35 GB of `Peptide`s. 25 s total.
- **Unspecific human, 7–25 aa:** `database_size.md`'s immunopeptidomics config, proteins
  checked on every 1000th row. Every one of the 673,210,816 rows materializes back
  exactly; 0 of 673,222 protein checks differ.

  | range | rows | compact | monolithic index + `Peptide`s (`database_size.md`) |
  |---|---:|---:|---:|
  | class I 8–12 | 148.9M | 2.38 GB | 42.4 GB |
  | class II 9–25 | 620.3M | 9.93 GB | ~275 GB |

  The build's peak RSS was 31.9 GB at length 25: sage's own single-length digest
  plus the runs accumulated so far. The merge peaked at 38.2 GB (the runs plus the
  merged table) and found 0 rows out of order. That run used the single-threaded
  map and longest-run pivots. The final code produces the same smallest offsets,
  and was re-checked:
  - unspecific 7–12: 174,273,871 rows, 0 mismatches, merged in 1.96 s, 0 out of
    order;
  - tryptic `--full`: unchanged.
