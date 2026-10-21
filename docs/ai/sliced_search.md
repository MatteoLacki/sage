# Sliced search: one precursor-mass slice of the index at a time (2026-10-13)

Stage 2 of `plans/sliced_search.md` (necromerge2 root). With `database.sliced: true`,
`Runner::new` builds a `CompactPeptides` table (`compact_peptides.md`) instead of an
`IndexedDatabase`. `sliced_search::search` then builds and searches the fragment index
one slice of that table at a time. An unspecific human class I digest needs 42 GB as a
monolithic index; sliced, it fits in a few GB (`database_size.md`).

Since 2026-10-21 the sliced search is the default (`Input::build`): `database.sliced`
unset means sliced unless one of three monolithic-only options is set, and `sliced: true`
with any of them is rejected:
- **`prefilter`:** another way of bounding memory; they don't compose, and sliced is the
  exact one (prefilter's first pass ranks candidates within protein chunks).
- **`wide_window`:** each slice is extended by the widest candidate window, which
  `wide_window` makes as wide as the isolation window, and every charge of a spectrum is
  pooled into one candidate list, which sliced splits per charge.
- **`chimera`:** each further PSM is found after removing the previous PSM's peaks, so
  the result depends on search order; sliced runs that loop per charge. Not supported on
  purpose: deconvolving co-isolated peptides (NNLS or similar) is the better route.

Measured on the F9477 best job (tryptic, one slice, same binary, 2026-10-21; necromerge2
`temp/peptide_views/sliced_best/`): peak anonymous memory 6.75 GB monolithic, 3.84 GB
sliced; wall 149.4 s monolithic pages, 161.4 s sliced pages, 106.4 s sliced grid (the
grid makes the speed, slicing the memory). SAGE 1% spectra / peptidoforms / ions
95,216 / 21,046 / 23,031 monolithic, 95,262 / 21,135 / 22,975 sliced; mokapot 122,144 /
28,884 / 31,680 and 119,471 / 28,661 / 31,343 (per-charge trims, see below).

## Why a slice gives the same answer as the whole index

A query (spectrum, precursor charge) only ever reads peptides whose mass falls in
its candidate window: the hull, over isotope errors, of
`effective_precursor_tol.theoretical_range(M − e·NEUTRON)` (`Scorer::
candidate_mass_window`, built from the same `searched_isotope_errors` and
`candidate_charges` that `initial_hits` and `matched_peaks` use). Rows are sorted by
mass, so the window is a contiguous row range.

Exactness needs more than every candidate being present. `IndexedDatabase::query`'s
`binary_search_slice` returns `[last row below the window, first row above it)`, one
padding row included. `matched_peaks_with_isotope` sizes its preliminary array from
that range, and `bounded_min_heapify` trims it. The trimmed order decides
`build_features`' stable-sort ties: equal rank keys, e.g. I/L isomers, keep
preliminary order. So a slice must reproduce that exact range, shifted. `SlicePlan`
(`slicing.rs`) guarantees it:
- **Cores split only where the mass strictly increases.** A query goes to the core
  whose first mass is at or below its window's low end. So the global "last row below
  the window" is at or after the row before the core.
- **Each slice holds:**
  - that one padding row below the core;
  - the core;
  - every row lighter than the next core's first mass plus the widest window (plus
    0.01 Da against float rounding).

  The window's high end is below that limit, so every row the global search reaches
  is in the slice.

`slicing::test::every_window_finds_its_global_rows_inside_its_slice` checks this on
random windows and on windows starting exactly on core boundaries: the local
`binary_search_slice` equals the global one, shifted.

Since Stage 3 of `plans/streaming_spectra.md` (necromerge2 root) a slice is built over
fewer rows: `slicing::rows_reaching` takes, for the queries actually assigned to it,
from one row below their lowest window start to the first row above their highest
window end (both widened by the 0.01 Da margin). Any window of those queries finds the
same `binary_search_slice` range in those rows as in the whole table, padding row
included, and they lie inside the planned slice above
(`windows_find_their_global_rows_inside_the_rows_their_slice_reaches`). The cores, and
so which queries share a slice, are unchanged.

Inside a slice, peptide `i` is table row `first_row + i`. Several things hold under
that shift:
- **`PreScore` order is unchanged:** it orders by matched count first, then
  `PeptideIx`, and the shift keeps `PeptideIx` order. The `u32::MAX` sentinel only
  sits on entries with 0 matches.
- **Hit counting is order-independent:** integer counts per candidate.
- **The page layout doesn't matter:** a different layout changes how pairs are
  found, not which (the same reason the `bucket_size` change was output-exact,
  `simd.md`).
- **Positional predictions line up:** they are passed as subslices at `first_row`
  (`Runner::make_scorer`, scaled by the charge span for IIM and fragment intensity).

## Charges

`plan_queries` makes one query per charge `initial_hits` would try
(`Scorer::candidate_charges`): every charge listed for the first precursor's ion,
including multi-charge pmsms precursors (charges 23, 24, 34, 234: 26% of F9477's
precursors), which until 2026-10-15 were searched at their first charge only (see
`multi_charge_precursors.md`).

A spectrum searched at one charge is fully exact. A spectrum searched at several
(listed charges, unannotated, or `override_precursor_charge`) becomes one query per
charge, scored on a copy whose only precursor carries that charge (slice scorers run
with `override_precursor_charge = false`). Those queries can land in different slices, so
`rank_across_charges` merges them:
- **Re-ranked together:** the best `report_psms` by rank key, with `rank` and
  `delta_best` recomputed. `delta_next` is taken against the next reported PSM or
  each charge's best unreported candidate (its last PSM's key minus that PSM's
  `delta_next`), whichever is higher.
- **Per charge instead of pooled:** `poisson`'s λ, `scored_candidates`, and the
  preliminary top-k trim. That is the accepted difference.

A pooled search trims to one top k across all charges; per-charge queries keep k
each, so they fully score a superset. Measured on G23 with
`override_precursor_charge: true` (every spectrum searched at charges 2–4), sliced
against monolithic:
- **Same winner:** 21,968 of 23,181 spectra, but `poisson`, `scored_candidates` and
  the LDA outputs differ on all of them.
- **Different winner:** 1,213 spectra. The sliced winner's hyperscore was never
  lower: higher in 1,045, tied in 168.
- **Different isotope error, same winner:** 6,629 spectra, all with identical
  hyperscores (isotope errors of one peptide tie on fragments).
- **`delta_next`:** smaller in 637 cases, because a better runner-up was scored.

## After the search

Features come back in spectrum order (rank order within a spectrum), as the
monolithic `par_iter` produces them. `reported_database` renumbers `peptide_idx`
from table rows to an `IndexedDatabase` holding only the reported peptides,
materialized with proteins and in row order. The rest of `run()` (LDA, FDR, protein
grouping, RT/IIM models, LFQ, writers) then runs unchanged. None of it reads beyond
reported peptides or writes raw `PeptideIx`.

## Slice size

A core's cost is about 200 B per row (its materialized `Peptide`) plus 16 B per
fragment (the index entry plus the sort and unzip copies while building). Cores are
filled up to half of `MemAvailable`, measured when the search starts.
`SAGE_SLICE_FRAGMENTS=<n>` caps fragments per core instead, for tests. Slices with no
queries are never built.

Since Stage 3 of `plans/streaming_spectra.md`, when every input can be read by index
(`indexed_spectra.rs`: pmsms, TDF, local mzML and MGF; not with TMT or LFQ) all files
are planned together and each slice is built once per run. Planning reads only each
spectrum's header (MS level, first precursor; pmsms from its precursor table); the
`min_peaks` check moved from planning to scoring, where a spectrum below it yields no
features either way. A slice reads the spectra of its queries when it is searched, one
read per spectrum even when several of its charges land there. Otherwise files are
read whole in SAGE's batches, and slices are rebuilt per batch. Each slice logs its
build and search time with the share of the cores each kept busy.

## Verification

Compared PSM by PSM on (filename, scannr, peptide, charge, label, rank), every
column except `psm_id` and `sage_discriminant_score`, plus all matched fragments:
- **F9477 production `run_sage` command** (predicted RT and fragment intensities):
  - new binary, monolithic: identical to the production node (361,854 PSMs, 2,282,445
    fragments);
  - sliced, `SAGE_SLICE_FRAGMENTS=20000000` (10 slices): identical to monolithic.
    Search 34.3 s against 37.8 s; total 67 s against 81 s;
  - sliced, automatic budget (1 slice): identical.
- **G23 DDA `.d`, tryptic human** (all precursors have a charge), 10 slices:
  identical (16,661 PSMs, 79,866 fragments). Reading `.d` first needed
  `util::local_path`: inputs are `file://` URLs, which timsrust can't open.

- **Unknown-charge merge path:** see Charges above.

## Scale (F9477 pmsms, unspecific human, automatic budget)

Config: F9477's search settings (±6.7 ppm per-precursor windows, isotope errors 0–3),
database `cleave_at = ""`, variable M-ox, protein N-term acetyl, cysteinylation,
`max_variable_mods = 1`, no positional predictions. F9477 is not an HLA sample, so
these runs measure memory and time, not identifications.

| | class I, 8–12 aa | class II, 9–25 aa |
|---|---|---|
| compact table | 148.9M rows, 2.38 GB, built in 144 s | 620.3M rows, 9.92 GB, built in 564 s |
| slices | 3 (19.6M–70.1M rows, 343M–950M fragments) | 22 (8.6M–47.8M rows, 376M–967M fragments) |
| slice build / search | 50 s / 140 s | 390 s / 348 s |
| whole run | 6 min 15 s | 22 min 39 s |
| peak RSS | 42.0 GB | 43.8 GB |
| monolithic index (`database_size.md`) | 42 GB resident, about twice that to build | about 275 GB |

In class II, slices 12–22 hold the heavy peptides only a few thousand F9477
precursors reach, but each still costs a full build (about 20 s, 220 s in all).
Narrowing a slice's rows to the windows its queries actually span would cut that;
not done.

The 42 GB peak follows from the budget: the slices are sized to half the memory
available when slicing starts (about 30 GB), and glibc keeps the compact build's freed
heap mapped. Fewer, larger slices search faster; `SAGE_SLICE_FRAGMENTS` trades speed
for a lower peak.

Since 2026-10-15 a non-specific table is generated directly (`compact_peptides.md`).
F9477 class I, same binary otherwise (`70017ca` against the direct build):

| | per-length build | direct build |
|---|---:|---:|
| compact table | 144.0 s | 7.7 s |
| anonymous memory when slicing starts | 11.8 GB | 2.4 GB |
| slices (fragments in the largest) | 3 (926M) | 3 (1,100M) |
| sliced search | 233.1 s | 234.5 s |
| whole run | 394 s | 259 s |
| peak anonymous / RSS | 35.4 / 46.3 GB | 41.5 / 52.5 GB |

Output identical (484,051 PSMs, 3,308,239 fragments). The higher peak is the budget:
no freed digest heap is left mapped, so more memory is available when slicing starts
and the slices come out larger; search time did not change.
