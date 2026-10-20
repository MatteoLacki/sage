# Merged fragment windows (`merge_fragment_windows`, 2026-10-20)

The preliminary score counts one match per (peak window, fragment) pair. Two peaks closer than
the fragment tolerance therefore both count a fragment lying between them, so a candidate's
prescore grows with near-duplicate peaks rather than with the fragments it explains. With
`merge_fragment_windows: true` (config, default `false`) the windows of a spectrum are merged
before the index is searched (`merge_overlapping_windows`, `scoring.rs`): each run of
overlapping or touching windows becomes one window over their union. A fragment then counts once
however many peaks lie within tolerance of it; a peak within tolerance of two fragments still
counts for both. Only the prescore's windows change; full scoring, which takes the most intense
peak per fragment, is untouched. The Poisson feature's lambda (prescore matches per candidate)
changes with it.

- `matched_peaks` builds a spectrum's windows once per precursor charge (`fragment_windows`) and
  hands the same windows to every isotope error's search (`matched_peaks_with_isotope`); they were
  rebuilt per isotope error before, which the merge would have repeated four times per query.
- Windows arrive in order of their lower bounds (peaks are in mass order, one window each) and
  are merged in one pass, in place, each window's bounds computed once; other orders (several
  fragment charges per peak) are sorted first.
- A merged window is stored as `(hi, Da(0, hi - lo))`. Its `theoretical_range` is exactly
  `(lo, hi)`: for bounds within a factor of two `hi - lo` is exact (Sterbenz), and then so is
  `hi - (hi - lo)`. Windows that overlap nothing are kept as they were.
- `merged_windows_hold_each_mz_once`: every m/z, window bounds and their neighbouring floats
  included, lies in one merged window when it lay in any window, and in no more than one.

## Measured (F9477, sliced grid, 1% FDR from SAGE's own q-values, no mokapot)

Prototype (necromerge2 `temp/peptide_views/merge_windows.patch`): about 9-10% of peak windows
overlap a neighbour (class I 1,271.5M -> 1,142.8M windows in its first slice, tryptic 2,701M ->
2,457M) and prescore matches drop 4.4-5.4%. Counts change in 74-96% of searches with a
candidate, and the top-50 set in 55-95% (class I's many equal counts reshuffle on ties). IDs
(spectra / peptidoforms / ions): class I 13,191 / 3,227 / 3,402 -> 13,312 / 3,245 / 3,417,
tryptic 78,748 / 19,327 / 21,114 -> 78,817 / 19,347 / 21,220.

Final SAGE of `jobs/f9477_best.toml` (inputs of `nodes_f9477_best_charges_20261015`, predictions
converted to mmappet, pages layout) and its mokapot step, same binary (necromerge2
`temp/peptide_views/merge/best/`):

| | spectra | peptidoforms | ions |
|---|---:|---:|---:|
| SAGE, off (= recorded run 17) | 94,053 | 20,705 | 22,717 |
| SAGE, on | 95,216 | 21,046 | 23,031 |
| mokapot, off (= recorded run 18) | 121,372 | 29,120 | 31,778 |
| mokapot, on | 122,144 | 28,884 | 31,680 |

SAGE's own counts rise 1.2-1.6%; after mokapot the change (-0.8% peptidoforms, -0.3% ions) is
within mokapot's run-to-run noise (1-2%).

Speed is unchanged (warm, alternating): best job search 115.9 / 112.0 s off, 112.6 / 114.6 s on;
class I grid 60.9 s (876 CPU-s) off, 59.0 s (889) on; tryptic grid 33.2 s (499) off, 34.0 s (506)
on. Off is identical to the binary before the flag in every check; on is identical to the
prototype.
