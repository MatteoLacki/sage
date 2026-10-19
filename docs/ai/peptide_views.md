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

## Stage 2: indexes and scoring over compact rows

`IndexedDatabase` has a row source: `peptides` (built `Peptide`s, as before) or
`compact: Some(CompactSlice)` (an `Arc<CompactPeptides>` and the slice's row range;
`peptides` then stays empty). `row_count`, `row_mass` and `peptide(idx) -> PeptideRef`
(`Built(&Peptide)` or `Compact(CompactView)`, itself a `PeptideView`) read either.
- `query` binary-searches `row_mass` (`binary_search_rows`, the same two partition
  points as `binary_search_slice`); `candidate_rows` reads `row_mass`.
- `Parameters::build_from_compact(table, rows, rt)`: the grid's `fragments_of` and
  the pages' fragment pass run `indexed_fragments` (now generic) over
  `CompactSlice::view(row)`; `index_pages`/`index_grid` are the shared tails of both
  builds.
- `Scorer`: `build_features` reads the `PeptideRef`'s few fields through the trait;
  `score_candidate` and `remove_matched_peaks` match the `PeptideRef` once and call
  generic bodies (`score_peptide`, `remove_peaks_matching`), so the per-residue loops
  are monomorphised for `Peptide` and for `CompactView`, not matched per residue.
- The sliced search builds every slice with `build_from_compact` and never calls
  `materialize_rows`; `Runner.compact` is an `Arc`. Reporting still materialises the
  reported rows (`reported_database`).

The monolithic search, prefilter and the reported database keep `Vec<Peptide>` and
the same code paths (generic functions instantiated for `Peptide`).

`crates/sage/tests/integration.rs` had not compiled since the predictions became
`PredictedValues` (`b76d235`, which ran only `--lib` tests); its RT/IIM tests now use
`f64` rows with NaN where they meant "no prediction" (every row has a value in
positional files; NaN never fits a window, so those candidates are evicted where the
tests' `None` kept them, and the tests' assertions are unchanged).

### Checks

- `compact::test::index_over_compact_rows_equals_index_over_materialized_rows`: pages
  and grid over rows 3..n-2 of the colliding-mods table, built from compact rows and
  from the same rows materialised: equal row masses, fragments, and `query` /
  `candidate_rows` for every row's mass at two tolerances.
- Suites: sage-core (lib 186, integration 17), sage-cloudpath, sage-cli.
- Data, against the pages baselines (necromerge2 `temp/peptide_views/check.sh`), all
  identical: small pmsms and G23 sliced, grid and pages, table built and stored;
  F9477 class I grid and pages, class II grid, stored tables.

### Measured (F9477, 16 threads, grid layout, stored table, warm)

The slice budget (half of `MemAvailable`) charged 200 B per row for the `Peptide`
the slice used to hold; rows now cost nothing, so `BYTES_PER_ROW` is gone and a
slice is sized by its fragments alone (16 B each, the pages build's cost; the grid
needs ~4.5 B).

| | before (stage 1) | views, old 200 B/row budget | views, final budget |
|---|---:|---:|---:|
| class I slices | 2 | 2 | 2 |
| class I peak anonymous | 20.46 GB | 5.79 GB | 9.23 GB |
| class I wall | 102.1 s | 87.1 s | 89.0 s |
| class II slices | 13 | 13 | 10 |
| class II peak anonymous | 20.48 GB | 7.81 GB | 9.76 GB |
| class II builds / searches | 116.5 / 166.2 s | 114.2 / 175.5 s | 100.7 / 160.3 s |
| class II wall | 348.4 s | 326.7 s | 297.0 s |

- Anonymous memory after a slice's search is now ~1.4-1.7 GB (spectra state and
  features), against 10-11 GB before.
- With the old slices the class II searches were ~5% slower (full scoring reads
  modifications through the view: var-site and static-rule lookups per residue);
  the final budget's fewer slices more than make up for it. Precomputing a view's
  modifications once per scored candidate is the obvious step if that cost matters.
- Class I pages layout (stored table): identical, peak anonymous 24.3 GB (39.1 GB
  for the pages run of `062a871`), as the pages build still holds 16 B per fragment.
- F9477 best job's monolithic final search (unchanged path): identical, warm search
  118.9 and 118.4 s against 120.3 s.
