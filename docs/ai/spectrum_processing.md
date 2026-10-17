# Spectrum processing cost and the top-N selection (2026-10-16)

`SpectrumProcessor::process` (`crates/sage/src/spectrum.rs`) turns each read spectrum
into the peaks the search uses. With the grid search (`fragment_grid.md`) it became a
large share of search time. Measured with temporary in-source cycle counters
(necromerge2 `temp/fragment_grid/measure/`), deisotoping on, `max_peaks` 800:

| | F9477 class I sliced | F9477 run 11 (monolithic tryptic) |
|---|---:|---:|
| spectra read | 0.82M (some twice: multi-charge queries in different slices) | 1.25M |
| raw peaks per spectrum | ~2,370 | ~2,450 |
| read + process, CPU-s | 332 of 1,033 search | 531 of 1,624 search |
| `deisotope` | 125 | 199 |
| sort of every deisotoped peak by intensity | 110 | 199 |
| take the top 800 heads, sort by mass | 44 | 51 |
| reading (tof → m/z, recalibration) | 53 | 81 |

`process_ms2` sorted every deisotoped peak by (intensity descending, m/z ascending),
dropped isotope peaks (`envelope.is_some()`), took the first `take_top_n` and
stable-sorted them by mass. `most_intense_in_mass_order` gives the same peaks in the
same order without the full sort:
1. drop the isotope peaks;
2. `select_nth_unstable_by` the top `take_top_n` by (intensity descending, m/z
   ascending), linear time;
3. sort those by (mass, then intensity descending, then m/z), which reproduces the stable
   mass sort's order among equal masses.

The output is identical except where two distinct peaks tie on both intensity and m/z
and straddle the cut-off; there the old unstable sort and the selection may keep
different ones (both are arbitrary but deterministic). The unit test
`selection_of_most_intense_heads_equals_the_full_sort` checks the two against each other
on 200 random spectra at several `top_n`, with equal masses and intensity ties broken by
m/z; breaking the m/z tie-break, the mass-sort tie-break, the cut-off or the head filter
each fails it. On the data checks (small inputs in every format, G23, F9477 runs 10 and
11, class I, class II) the outputs are identical.

`deisotope` itself is a two-pointer scan that recomputes its ppm tolerance in the inner
loop; making it faster while keeping every float comparison identical was not tried.

Effect (same binaries otherwise, warm page cache): F9477 run 10 grid search 19–20 s →
13.0 s, run 11 grid 84.6 s → 75.9 s; class I sliced light slice 14–16 s → 10.4 s
(whole class I search 68.0 → 63.8 s with compact table v2 and the RT code path, which
do not change it).

