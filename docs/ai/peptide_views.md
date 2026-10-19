# Peptide views: fragments and scoring without a `Peptide` (2026-10-19)

Plan: necromerge2 `plans/peptide_views.md`. A sliced search turns every row of a
slice into a sage `Peptide` (`CompactPeptides::materialize_rows`), ~200 B each with
its heap allocations: on F9477 the first slice's peptides were 16.1 GB (class I,
82.3M rows) and 15.7 GB (class II, 76.8M rows) of the slice's ~20.9 GB anonymous
memory, the grid and directory 4.4 and 4.7 GB. They are read only to generate
fragments while the index is built and, during the search, by full scoring of at
most 50 candidates per query (`trim_hits`).

## Stage 1: `PeptideView` and `CompactView`

`peptide_view.rs`: `PeptideView` is what fragment generation and scoring read of a
peptide: `len`, `residue(i)`, `modification(i)`, `nterm`, `cterm`, `monoisotopic`,
`decoy`, `missed_cleavages` (`label` derived). `impl PeptideView for Peptide` returns
the fields. `IonSeries<'p, P: PeptideView + ?Sized = Peptide>` is generic; every
existing call passes a `&Peptide` and is monomorphised to the same code as before.

`compact/view.rs`: `CompactPeptides::view(row)` is a `CompactView`, the decoded
16 B row, the table and the target's residue window in the proteome copy; nothing is
allocated. It reads as `materialize` builds, which is sage's own code, so it follows
sage's rules rather than an obvious "sum of mods":
- **First mod wins.** `Peptide::apply_site` puts a mass on a residue only while it
  holds 0, on a terminus only while it is `None`; `static_mods`' residue rule assigns
  only where 0. `materialize` applies the variable sites first, then the static mods
  in the iteration order of the run's `static_mods` `HashMap`. So a residue with a
  variable mod gets no static mod, and when static rules compete for a residue (e.g.
  `^M`, `[M` and `M` on a protein-N-terminal M) the first in that order wins.
  `StaticMods` captures that order once, from the same `HashMap` instance the table
  holds (built in `CompactPeptides::new`), plus a per-residue table for interior
  residues, where only residue rules apply.
- **Protein-terminal rules** (`[X`, `]X`, `[`, `]`) apply only when the row's
  position is protein N/C-terminal or the whole protein.
- **Decoys** are computed on the target and read through `Peptide::reverse`'s index
  map: first and last residues kept, interior `i -> len - 1 - i` (for 3+ residues).
- `apply_site` adds to 0 (`0.0 + mass`), a residue rule assigns `mass`; the view does
  the same, and a rule leaving 0 does not block later ones.
- The mass is the row's stored `monoisotopic` (bit-identical to sage's,
  `verify_compact_peptides`).

### Checks

- `compact::test::views_equal_materialized_peptides`: every row of the compact tests'
  digests (tryptic, semi-tryptic, non-specific with terminal mods, the edge-case FASTA
  with every non-specific config), plus two configs made for this: static rules
  competing for the same residue and terminus with variable sites on the same
  residues, and protein-terminal residue rules with nothing competing.
  `peptide_view::first_difference` compares residues, modification bits, both
  termini, mass, decoy, missed cleavages and every ion of all six kinds bit for bit.
  Applying statics in reverse order, summing a variable and a static mod, dropping
  the decoy reversal, and ignoring the protein position (N residue, C residue, N
  terminus) each fail it; the protein-position mutations were caught only after the
  non-competing config was added (with competition, another rule usually wins first).
- `verify_compact_peptides --views N` (F9477, human): production tryptic config
  (static C, variable M and protein N-term acetyl) all 6,300,882 rows; class I every
  2nd row (74.5M of 148.9M); class II every 8th (77.5M of 620.3M): 0 differences.
- Suites pass; the F9477 best job's final search with the generic `IonSeries` is
  identical to the previous binary's, and as fast within noise (warm 118.0 s vs
  116.6 s; the round before 120.7 vs 124.7 s).
