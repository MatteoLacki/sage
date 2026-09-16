# ms2_* mokapot feature regression — found and fixed (2026-09-15)

## Symptom

A fresh, code-correct rerun of the F9477 "best known" job (`jobs/f9477_best.toml`)
yielded ~20% fewer PSMs/peptides at 1% FDR than the stored 2026-09-07 baseline
(99,112/28,275 → ~79k/~22k in the broken run). MS1, clustering, precursor counts,
and raw hit counts were all byte-identical or near-identical between old and new —
the deficit traced specifically to mokapot's `ms2_*` fragment-intensity-similarity
features (entropy, cosine, pearson, etc.) collapsing toward zero/noise.

## Root cause

Commit `5753d07` ("Simplify score_candidate's Option-heavy fragment-intensity
plumbing", 2026-09-10) extracted `Scorer::iter_fragments` from inline code in
`crates/sage/src/scoring.rs`. In doing so it moved `.enumerate()` from *inside*
the `ion_kinds.iter().flat_map(...)` closure (per ion kind, `idx` resets to 0 at
the start of each kind's `IonSeries`) to *outside*, wrapping the whole chained
iterator in `score_candidate`:

```rust
// before (correct): idx resets per kind
.flat_map(|kind| IonSeries::new(peptide, *kind).enumerate())

// after (bug): idx runs globally across kinds
.flat_map(move |kind| IonSeries::new(peptide, *kind))
// ... later: self.iter_fragments(peptide).enumerate()
```

With `database.ion_kinds = ["b", "y"]` (the real config), every b-ion consumed
one unit of the shared counter before the y-series started, so every y-fragment's
`idx` was offset by `peptide.sequence.len() - 1` (the b-series length). This
corrupted two things:

1. The y-ordinal reported in `matched_fragments.sage.tsv`
   (`len - 1 - idx`), producing impossible negative ordinals (e.g. `y -11` for a
   14-residue peptide) — the first visible symptom that cracked the case.
2. `ms2_similarity::fragment_annotation_id(Kind::Y, idx, charge)`, the slot used
   to look up predicted MS2 intensity for comparison against the observed
   spectrum — silently comparing against the wrong (or out-of-range) predicted
   intensity for essentially every y-ion match, collapsing every `ms2_*` mokapot
   feature.

The commit's own claim ("no behavior change — verified against full sage-core
unit and integration suites") was wrong: no test exercised ordinal correctness
across multiple chained ion kinds.

## How it was found

1. Ruled out (with direct evidence, not inference): compact-MS1 migration,
   dump_peptides mass-tie sort ordering, an in-progress `PeakColumns`
   struct-of-arrays refactor (uncommitted WIP in this repo — proved innocent by
   `git stash`-ing it and reproducing byte-identical broken output with vs.
   without it), mz-sort cache corruption, fragment-intensity cache
   coverage/pointer-mapping correctness, charge-range consistency.
2. Localized the regression specifically to `ms2_*` features via mokapot PIN
   feature-column statistics (old baseline `.pin` vs. broken run's `.pin`).
2. Found a real PSM with `has_predictions=true`, a verified-correct predicted-
   intensity pointer, yet near-zero `ms2_cosine_similarity` — and noticed
   `matched_fragments.sage.tsv` reported negative y-ordinals for it.
3. Confirmed the old 2026-09-07 baseline's `matched_fragments.sage.tsv` has
   zero negative ordinals anywhere (2.3M rows checked) — so this was a genuine
   regression, not a pre-existing wart.
4. `git log` on `scoring.rs`/`ion_series.rs` since the baseline date turned up
   `5753d07`; its diff showed the `.enumerate()` relocation directly.

## Fix

Commit `548ffd5` (branch `optimizations/simd`, this fork). `iter_fragments` now
returns `impl Iterator<Item = (usize, Ion)>` with `.enumerate()` back inside the
per-kind closure; both call sites (`remove_matched_peaks`, `score_candidate`)
updated accordingly. This commit is isolated from the pre-existing uncommitted
`PeakColumns` WIP also present in this working tree (kept separate via
stash → fix on clean HEAD → commit → stash pop; verify the WIP diff doesn't
reintroduce the bug before touching `scoring.rs` again).

## Verification

- `cargo test --release -p sage-core --lib`: 143/143 (145/145 with the WIP
  layered back on).
- Same exact `sage` invocation (fasta/pmsms/precursors/predicted-rt/
  fragment-intensity-cache all held fixed), before vs. after fix, for the
  whole F9477 dataset:

  | | ms2_cosine mean | frac zero | ms2_pearson mean | ms2_entropy mean |
  |---|---|---|---|---|
  | old baseline (2026-09-07) | 0.4345 | 0.004 | 0.2781 | 0.4795 |
  | broken (pre-fix) | 0.1760 | 0.115 | 0.0306 | 0.2285 |
  | **fixed** | **0.4349** | **0.004** | **0.2788** | **0.4796** |

- Full necroflow pipeline rerun (`jobs/f9477_best.toml`) through mokapot:
  PSMs/peptides @ 1% FDR went from a ~20% deficit to 97,143/27,889 — within
  ~2% of the 99,112/28,275 baseline (residual gap consistent with normal
  mokapot/xgboost training stochasticity).

## Incidental notes for whoever picks this up

- necroflow's `run_sage` (and its dependents like `recalibrate_pmsms_mz`) can
  re-execute **in place** at the same node-hash directory when a `source_*`
  binary's content changes but its recorded `node_key`/config string doesn't —
  some downstream rules assume a clean/nonexistent output path and will raise
  `FileExistsError` on rerun (`searchops/recalibration.py`'s
  `recalibrate_pmsms_mz`). Had to manually delete four stale output files in
  `nodes_f9477_best_run2/recalibrate_pmsms_mz/<hash>/` before the pipeline
  would proceed. Not fixed here — just worked around for this verification run.
- The `PeakColumns` refactor in this working tree (`spectrum.rs`, and its
  ripple into `scoring.rs`/`runner.rs`/`lfq.rs`/`tmt.rs`/`tests/integration.rs`)
  is unrelated, pre-existing, uncommitted work — not authored as part of this
  investigation. Left untouched/uncommitted, exactly as found.
