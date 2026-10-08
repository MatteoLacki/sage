# Tolerance convention: windows on `observed − theoretical` (2026-10-08)

Every `Tolerance` (`crates/sage/src/mass.rs`) is a window `[lo, hi]` on
`observed − theoretical`, where "theoretical" is whatever the observation is
compared against: a peptide/fragment mass, or an externally predicted RT/IIM.
This covers `precursor_tol`, `fragment_tol`, `rt_tol_sec`, `mobility_tol` and the
pmsms per-precursor `ppm_tol_lo/hi`. It matches how the pipeline fits windows:
every residual it fits is `observed − theoretical` — RT (`rt − koina_rt`) and
IIM (`ion_mobility − koina_iim`) in `git/featureprediction`, precursor
(`precursor_ppm`, exp − calc) and fragment (`fragment_ppm`) in `git/searchops` —
and SAGE's own reported `delta_mass`/`precursor_ppm` uses the same sign.

## API

There is no `bounds(center)`. Each call site says which side it knows:

- `theoretical_range(observed)` = `[observed − hi, observed − lo]`: candidates
  compatible with an observation.
- `observed_range(theoretical)` = `[theoretical + lo, theoretical + hi]`:
  observations compatible with a candidate.

Ppm/Pct windows scale by whichever value is known. Against SAGE's symmetric-
denominator `delta_mass` the difference is O(ppm²) relative (~1e-10), below f32
resolution (~6e-8), so a PSM's reported `precursor_ppm` lies inside the window it
was searched with.

| Site | Known value | Method |
|---|---|---|
| `database.rs` `query` (precursor) | observed precursor mass (minus isotope error) | `theoretical_range` |
| `database.rs` `page_search`/`page_search_batch` (fragments + precursor filter) | observed peak mass | `theoretical_range` |
| `scoring.rs` `evict_rt_iim_mismatches` | observed RT / IIM | `theoretical_range` |
| `spectrum.rs` `select_matched_peaks` (final fragment matching, TMT) | theoretical fragment / reporter | `observed_range` |
| `spectrum.rs` `in_isolation_window` | target m/z | `theoretical_range` |
| `lfq.rs` MS1 mass / mobility ranges | calcmass / PSM ims | `observed_range` |

Isolation windows are physical ranges, not residual windows, so readers convert
them: mzML lower/upper offsets become `Da(-upper, lower)` (candidates in
`[target − lower, target + upper]`); tdf and mgf windows are symmetric.

## Why: the bug this replaced

Upstream `Tolerance::bounds(center)` returned `center + [lo, hi]` — correct only
when `center` is the theoretical value. The precursor query, the preliminary
fragment page search and the RT/IIM eviction all passed the *observed* value, so
they applied every fitted window mirrored (`theoretical − observed`), while final
fragment matching (`select_matched_peaks`, centred on the theoretical fragment)
used the fitted sign. The same fragment window was therefore read both ways round
by the two search stages. Every test used symmetric windows, so nothing noticed.

Impact before the fix was set by how far each window's centre sat from zero.
F9477 `correct_precursors_rt` binned RT windows had centres from −1.6 s to +4.2 s
against half-widths of 6–25 s (the XGB correction uses mz/iim/charge, not RT, so
per-RT-bin medians are not zero): mirroring shifted each window by ~0.7–1σ.
Precursor windows were centred ~0.07 ppm against ±6.7 ppm, so negligible.

## Divergence from upstream

Upstream configs with asymmetric windows mean the opposite here: upstream's
open-search `precursor_tol: {"da": [-500, 100]}` is `[-100, 500]` in this fork.
The upstream config warning "Tolerances are applied to experimental masses, not
theoretical" encoded the old convention and was removed. Every hand-written window
in `necromerge2/jobs/*.toml` was symmetric at the time of the change.

## Verification on F9477 (`jobs/f9477_best.toml`, 2026-10-08)

The final search had byte-identical inputs before/after (config, RT-corrected
precursors, predicted RT); only the binary differed. Replicating SAGE's f32
arithmetic per PSM:

| Final search | outside mirrored window | outside intended `obs − theo` window |
|---|---|---|
| precursor, old binary | 0 | 3,196 |
| precursor, new binary | 2,754 | 0 |
| RT (`(rt − predicted_rt_external)·60` vs `rt_tol_sec(rt)`), old | — | 11,634 |
| RT, new | — | 0 |

Unique target peptides / ions at 1% mokapot FDR: 28,583 / 34,995 → 28,612 / 34,928
(+0.10% / −0.19%, within mokapot's input-sensitivity noise). The first,
symmetric-window search gives the same 207,506 PSMs with all values equal; only
`psm_id` (and sometimes `sage_discriminant_score`) differ, and they differ just the
same between two runs of the *old* binary, since neither is run-to-run deterministic.

Runtime, first search (identical work), 6 ABBA-ordered runs per binary: wall
36.19 ± 0.47 s → 36.44 ± 0.50 s (not significant); search stage 12.36 → 12.51 s;
"file IO" stage 5.56 → 5.75 s, consistently. The IO path's source is unchanged, so
the IO shift is attributed to fat-LTO code-layout changes. That is unverified:
instruction counts need `perf`, and `perf_event_paranoid` was 4. Final search
search-stage: 35.9/37.0 s → 36.6/36.1 s.

## Tests

`mass.rs`: exact asymmetric Da/Ppm/Pct ranges, and `t ∈ theoretical_range(o)` ⇔
`o ∈ observed_range(t)`. `tests/integration.rs`: asymmetric RT, IIM, precursor and
fragment windows, each keeping a candidate inside the window and rejecting one at
the mirrored offset. The fragment test needs both stages to agree. With
`theoretical_range` mutated back to the old `center + [lo, hi]`, all four
integration tests and both `mass.rs` tests fail.
