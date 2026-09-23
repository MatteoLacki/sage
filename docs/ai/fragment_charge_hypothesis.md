# Fragment charge-hypothesis loop vs. resolved peak charge

`matched_peaks_with_isotope`, `remove_matched_peaks`, and `score_candidate`
all loop a candidate fragment charge over `1..max_fragment_charge`, rescaling
one side of the comparison per iteration (see the "Experimental peaks are
multiplied by charge, therefore theoretical are divided" convention: whichever
array is the fixed, repeatedly-searched one keeps its native scale; the
per-iteration side gets rescaled). Algebraically `mass * charge == M_frag` and
`mass == M_frag / charge` are the same equation, so which side gets scaled is
a performance choice, not a correctness one -- these two forms are equivalent,
full stop.

## What is/was not a bug

An earlier read of this code (2026-09-23 chat) mistook the multiply-vs-divide
difference between `matched_peaks_with_isotope` (peak side) and the other two
functions (fragment side) for a discrepancy -- "one got fixed, one didn't".
That's wrong: `7cd21cd` (the commit that reworded the comment in the other two
functions) states "No behavior change" in its own message, and the algebra
above confirms it. Retracted; noted here so it doesn't get re-litigated.

## What is a real imprecision

`peak.mass` (`spectrum.rs`) is built as `(mz - PROTON) * scale_charge`, where
`scale_charge` is the peak's own resolved charge when deisotoping
(`deisotope: true`, the default) found an isotope envelope for it, else `1`
as the best available base scale when it didn't. When the charge really is
resolved, `peak.mass` **already is** that peak's true, charge-independent
neutral mass. There is no remaining charge hypothesis to test for it: the
loop's `charge` variable is answering "what if I apply *another* charge on
top of this peak's own, already-resolved one", which is not a real physical
question once the charge is known. Only the unscaled comparison (`charge ==
1` in the loop, i.e. no further rescale) is meaningful for such a peak; every
other `charge` iteration tests a fragment mass at a physically baseless
multiple/divisor of the peak's real mass, and can only match another
candidate's fragment by coincidence.

Consequence: those extra iterations can inflate a wrong candidate's
preliminary matched-peak count (in `matched_peaks_with_isotope`), which feeds
`quick_score`'s `min_matched_peaks` gate and `trim_hits`'s top-K truncation
*before* real hyperscore ever runs -- so in principle a spurious match could
crowd a genuinely correct candidate out of the trimmed set. Final hyperscore
(`score_candidate`) itself is unaffected for whichever candidates do survive
trimming; the same imprecision there and in `remove_matched_peaks` only risks
spurious matched-fragment bookkeeping for candidates already in hand, not
scoring correctness.

This is not fork-specific: verified identical in a fresh clone of upstream
`lazear/sage` (`master`, commit `2c9922e`, cloned 2026-09-23) -- same
peak-mass construction (`spectrum.rs`), same unguarded charge loop in all
three call sites. Upstream does not appear to treat it as a defect.

## Round 1 fix and its gap (2026-09-23, superseded below)

The first pass rejected `charge != 1` matches only when the peak's resolved
charge was `>= 2`, because `peak_charges[i] == 1` was, at that point,
irreducibly ambiguous: `deisotope()`'s raw output (`Deisotoped.charge:
Option<u8>`) does distinguish "confirmed charge 1" (`Some(1)`, envelope found
with `NEUTRON/1` spacing) from "no envelope found" (`None`), but
`process_ms2` collapsed both to the stored value `1` via `.unwrap_or(1)`
before it ever reached `ProcessedSpectrum`. So a genuinely-confirmed
charge-1 peak still got the full, imprecise `1..max_fragment_charge` loop --
the round-1 fix only closed the `>= 2` half of the problem.

## The real fix: a `0` sentinel for "unresolved" (2026-09-23, round 2)

`peak_charges[i]` now means "the resolved charge, or `0` if none was found" --
not "the charge, defaulting to 1". `process_ms2`'s deisotope branch
(`spectrum.rs`) keeps using `charge.unwrap_or(1)` to build `peak.mass`'s
*scale* (that part was never wrong -- 1 is still the right base to multiply
candidate charges onto when the charge is genuinely unknown), but now reports
`charge.unwrap_or(0)` as the stored charge instead of reusing that same
`unwrap_or(1)` value. This makes "confirmed charge 1" (`Some(1)`, stored as
`1`) and "no envelope, unknown" (`None`, stored as `0`) distinguishable for
the first time, closing the round-1 gap:

- `matched_peaks_with_isotope`: per-peak window construction now emits the
  unscaled window (`mass`, not `mass * charge`) for `peak_charges[i] >= 1`
  (any resolved charge, including confirmed 1), skipping the rest of the
  `1..max_fragment_charge` range entirely. Falls through to the full loop
  only for `Some(0)` or `None` (deisotope off, `peak_charges` empty).
- `remove_matched_peaks` / `score_candidate`: same "can't know which peak a
  search lands on in advance" reasoning as round 1, so the found peak index
  is still validated after the fact --
  `.filter(|&i| charge == 1 || query.peak_charges.get(i).copied().unwrap_or(0) == 0)`
  -- reject whenever `charge != 1` and the peak's stored charge is nonzero
  (resolved). `score_candidate`'s diagnostic `peak_charge` read (used for
  `exp_mz`/`calc_mz`/`fragments_details.charges`) needed its own fix: it must
  recover the *scale* charge peak.mass was actually built with (1 for
  unresolved), not the raw `0` sentinel -- `Some(0) | None => 1, Some(z) =>
  z` -- dividing by a literal `0` there would have produced `inf`/`NaN`.

Effect on this monorepo's actual jobs: none, same as round 1. Every
`jobs/*.toml` here pins `max_fragment_charge = 1`, so the loop is always
`1..2` -- a single iteration at `charge == 1` -- regardless of any peak's
resolved charge, before and after both fixes. This only matters if/when a job
sets `max_fragment_charge > 1` with `deisotope: true` (the default).

Verified: `cargo test -p sage-core --lib` (151 passed, up from 148 baseline --
three targeted tests, one replacing a round-1 test whose fixture semantics
changed under the sentinel), `cargo test -p sage-core --test integration` (12
passed), `cargo test -p sage-cli` (2 passed) -- all pre-existing tests
unchanged, no regression. Tests in `scoring.rs`'s `mod tests`:

- `resolved_multiply_charged_peak_only_tests_charge_one`: `peak_charges:
  vec![2]` at mass `b1/2` no longer spuriously matches `b1` at the
  `charge==2` iteration (`score.matched_b == 0`). Unchanged by round 2 (charge
  `2` meant "resolved" under both rounds).
- `confirmed_charge_one_peak_only_tests_charge_one`: same setup with
  `peak_charges: vec![1]` (now unambiguously "confirmed charge 1") -- also
  `score.matched_b == 0`. This is the exact case round 1 couldn't fix, and
  the reason round 2 exists. Confirmed it actually exercises the fix (not
  just compiles against it) by temporarily reverting the filter back to the
  round-1 `unwrap_or(1) < 2` predicate and re-running it in isolation: it
  failed with `matched_b == 1` as predicted; the revert was then undone and
  the full suite (151/151) re-verified green.
- `unresolved_peak_still_tests_multiple_fragment_charges`: same setup with
  `peak_charges: vec![0]` (genuinely unresolved, replacing round 1's
  `vec![1]` fixture, which stopped meaning "unresolved" under the sentinel)
  still matches at `charge==2` (`score.matched_b == 1`) -- confirms the fix
  doesn't over-reject the legitimate charge-hypothesis case.

All three only cover `score_candidate`'s filter (the same predicate is used
verbatim in `remove_matched_peaks`). `matched_peaks_with_isotope`'s window-
construction fix is simpler (skip building extra windows vs. post-hoc
filtering a found match) and mathematically the same claim, but isn't
independently covered by a test that reaches `self.db.query(...)` -- would
need a full precursor-tolerance-matching `IndexedDatabase` fixture; left
uncovered by a dedicated test here as a reasonable stopping point, not
because it's less important.
