# Closest-peak columns in the matched-fragments export

`matched_fragments.sage.tsv`/`.parquet` (written when `--annotate-matches` is set) carry two
extra raw m/z columns per matched fragment, alongside the pre-existing
`fragment_mz_calculated`/`fragment_mz_experimental` (which is always Sage's actual match — the
*most intense* peak in the tolerance window, unchanged): `closest_fragment_mz_calculated` and
`closest_fragment_mz_experimental`, the theoretical/observed m/z of the peak *closest in mass*
to the theoretical fragment within that same window, which may be a different, less intense
peak. When they're the same peak — the common case — all four columns line up
(`fragment_mz_experimental == closest_fragment_mz_experimental`); there is no null/
special-casing for that.

**This started as two signed-ppm-error columns instead** (`most_intense_fragment_ppm_error`/
`closest_fragment_ppm_error`), computed with a `symmetric_ppm_error` helper. That was dropped:
`most_intense_fragment_ppm_error` was a pure, deterministic function of two columns already
exported (`fragment_mz_calculated`/`fragment_mz_experimental`) — storing it was storing `A`,
`B`, and `A-B`. Raw m/z values are strictly more useful downstream too: a consumer can compute
ppm, Da, or anything else, in whatever sign convention it wants, instead of being locked into
one baked-in formula (and `searchops.recalibration` already computes its own ppm from raw m/z
columns — see below). If you're looking for `symmetric_ppm_error`/the ppm-error `Fragments`
fields, they no longer exist; this file's git history has the prior design if it's ever wanted
back.

## One scan, not two: `select_matched_peaks`/`PeakMatch`

`spectrum.rs`'s `select_most_intense_peak` became `select_matched_peaks`, returning
`Option<PeakMatch { most_intense: usize, closest: usize }>` instead of `Option<usize>`. Both
selections answer to the same tolerance window, so they're computed in one linear scan rather
than two near-identical functions. All 4 pre-existing call sites (`remove_matched_peaks`,
`observed_isotope_ladder`, `score_candidate`, `tmt::find_reporter_ions`) were updated; three of
them only use `.most_intense` and discard `.closest`.

Tie-breaks, both mirroring the original function's own `>=` "last write wins" policy:
- most-intense: unchanged — on an exact intensity tie, the higher-index peak wins.
- closest: smaller distance wins; on an exact distance tie, the more intense peak wins; on a
  tie of both, the higher-index peak wins.

## Charge-hypothesis filter applies to `closest` too

See `docs/ai/fragment_charge_hypothesis.md`: `score_candidate`'s per-charge loop only tests
`charge != 1` against peaks whose deisotoped charge is unresolved (`peak_charges[i] == 0`) — a
resolved-charge peak is only a legitimate match for the `charge == 1` hypothesis. The raw
`closest` index from `select_matched_peaks` is filtered by the identical rule before being
reported, falling back to the already-matched `most_intense` index if the geometrically closest
peak turns out to be charge-invalid for the hypothesis being tested. In practice this only
matters when `max_fragment_charge > 1`; every job in this monorepo pins `max_fragment_charge: 1`.

**Real-data confirmation (2026-09-23, `tests/config.json` smoke run, `deisotope: true`,
`max_fragment_charge: 1`)**: even at `max_fragment_charge: 1` — where the charge-hypothesis
filter is a no-op, since `charge == 1` always passes it regardless of a peak's resolved charge —
`closest` genuinely diverged from `most_intense` on 2 of 22 real matched fragments. In one case
the closest-by-neutral-mass peak was a peak deisotoped to a *different* resolved charge (2+)
than the matched peak, whose rescaled neutral mass coincidentally landed inside the same narrow
(±10ppm) window as the charge-1 theoretical fragment even though its own real observed m/z was
in a completely different region (~470 vs. ~938). This is not a bug introduced here — Sage's own
pre-existing most-intense matching already treats any resolved-charge peak as eligible at the
`charge == 1` test (that's what `docs/ai/fragment_charge_hypothesis.md`'s fix explicitly
permits) — but it means `closest_fragment_mz_experimental` can occasionally point at a peak in
an entirely different m/z neighborhood, not just "the next peak over." Downstream consumers
treating `closest` as "a small correction to `fragment_mz_experimental`" should be aware it can
be a large jump.

## Downstream consumer: `searchops.recalibration`

`git/searchops/src/searchops/recalibration.py` reads both m/z pairs by column name (see that
repo's `AI.md`, "`recalibration.py` — most-intense vs. closest fragment peak"):
recalibration-model *fitting* uses only fragments where `most_intense == closest`
(`_unambiguous_fragment_ppm`, since a coincidentally-close but physically different peak — see
the real-data case above — would bias the fit), while the *final reported ppm distribution*
(`plot_recalibrated_ppm`) uses `closest` unconditionally (`_closest_fragment_ppm`), since it
should reflect true mass accuracy rather than which peak happened to be most intense. Neither
function is defined here; this is just the contract the exported columns support.
