# Predicted MS2 fragment intensity: feature-only, no hard filtering (2026-08-31)

Third externally-predicted-property input, after `--predicted-rt`/
`--predicted-iim` (see `predicted_rt_iim.md`) — but a different shape by
deliberate design, decided explicitly when this was scoped: **feature-only,
no hard eviction filter**. `--predicted-rt`/`--predicted-iim` reject
candidates whose predicted value falls outside a tolerance window before
scoring finishes; fragment intensity never rejects anything. It only
computes one reported `Feature` field, `ms2_entropy_similarity`, from how
well a candidate's *predicted* fragment intensities (Prosit `compact_trt`,
via `git/featureprediction`) agree with its *observed* matched-peak
intensities — the kind of per-PSM MS2-similarity feature DIA rescoring tools
(DiaTracer, MSBooster) compute for downstream FDR refinement (mokapot/
Percolator), not something SAGE itself acts on directly.

## Two inputs, both optional, both-or-neither

- `--predicted-fragment-intensity-index <path.parquet>` — the small,
  job-scoped `(sequence, charge, start, end)` pointer index from
  `git/featureprediction`'s `export_fragment_intensity_for_sage` (see that
  repo's AI.md). `start`/`end` are a half-open row range into the second
  input.
- `--predicted-fragment-intensity-cache <path>` — the **shared,
  ever-growing** `arrays.mmappet` directory itself (a subdirectory of
  `git/featureprediction`'s fragment-intensity `PredictionCache` — this
  fork never touches that cache's `index.sqlite3`/`write.lock`, only its
  two flat arrays: `annotation_id` uint8, `predicted_intensity` float16).

`Input::build()` requires both or neither (`crates/sage-cli/src/input.rs`) —
unlike `predicted_rt`/`rt_tol_sec`/`rt_sigma_sec`'s three-way check, there's
no tolerance/sigma pair here at all, since nothing is evicted.

Real verified log line (2026-08-31, full F9477): `loaded 16874172 predicted
(sequence, charge) -> (start, end) fragment-intensity entries` /
`opened fragment-intensity arrays.mmappet (714530556 rows)`.

## Reader (`sage-cloudpath/src/fragment_intensity_cache.rs`)

- `read_fragment_intensity_index` — plain parquet row-iteration
  (`sequence`/`charge`/`start`/`end` columns) into a
  `HashMap<(String, u8), (u64, u64)>`. Same shape as
  `predicted_properties::read_predicted_iim`.
- `FragmentIntensityArrays::open` — mmaps `arrays.mmappet`'s two flat
  binary columns directly (`memmap2`, same technique `pmsms.rs` already
  uses for `pmsms.mmappet`/precursors `.mmappet`). Schema/mmap-slice
  helpers are a small, deliberate duplication of `pmsms.rs`'s private ones
  (different, incompatible error type; two call sites don't justify a
  shared abstraction). New `half = "2"` dependency (added to both
  `sage-cloudpath` and `sage-core`) reads the native FP16
  `predicted_intensity` column directly as `&[half::f16]` — already a
  transitively-resolved version in this workspace's `Cargo.lock` via
  arrow/parquet, so this added no new crate version to the dependency
  graph, just an explicit `Cargo.toml` line.

## Dense resolve: reuses `iim_dense_slot` unchanged

`Runner::resolve_predicted_fragment_intensity` (`sage-cli/src/runner.rs`)
mirrors `resolve_predicted_iim` exactly — one `Peptide::to_string()` call
per peptide, not per (candidate, spectrum) (same real ~37%-wall-time
regression this avoids, see `predicted_rt_iim.md`) — but reuses
`Scorer::iim_dense_slot` **as-is**, un-renamed: the `(peptide_idx, charge)
-> flat array slot` arithmetic is identical for any per-`(peptide, charge)`
dense array, not IIM-specific despite the name. Output:
`Vec<Option<(u64, u64)>>`, one slot per `(peptide_idx, charge)` in
`[min_precursor_charge, max_precursor_charge]`.

`Scorer` holds three new fields, all `Option`, all `None` when this feature
isn't configured (verified: `prefilter_peptides`'s `mini_runner` and both
`crates/*/tests/integration.rs` fixtures explicitly pass `None` for all
three — a lightweight prefilter pass never needs MS2 features regardless of
what the real search below it will use):

```rust
pub predicted_fragment_intensity_index: Option<&'db [Option<(u64, u64)>]>,
pub predicted_fragment_intensity_annotation_id: Option<&'db [u8]>,
pub predicted_fragment_intensity: Option<&'db [half::f16]>,
```

## Fragment-slot numbering: `ms2_similarity::fragment_annotation_id`

Must stay bit-for-bit consistent with `git/featureprediction`'s
`fragment_slots.py::annotation_id_for` — both sides were designed together
so the Rust reader needs **zero transform**:

```
kind_block(B) = 0, kind_block(Y) = 87   # 29 ordinals x 3 charges each
annotation_id = kind_block(kind) + idx * 3 + (charge - 1)
```

`idx` is SAGE's own raw backbone-position counter from
`IonSeries::enumerate()` — exactly what `score_candidate`'s existing
fragment loop already has (the *outer*, un-shadowed `idx`, before the
`annotate_matches`-only ordinal remap later in the same function) — no
ordinal recomputation needed on this path at all, unlike
`annotate_matches`'s own `fragment_ordinals` output. Fragments outside the
Prosit `compact_trt` vocabulary (kinds other than B/Y, `idx >= 29`, fragment
charge `> 3`) return `None` — same "no cache coverage" treatment as a
missing key, not an error.

## `score_candidate`: unpack once, accumulate in the existing loop, score once

Same three-step shape hyperscore itself already uses (accumulate
incrementally during one pass over the fragment ladder, compute the
derived score once after the loop — see the comparison drawn out loud
during implementation):

1. **Before the loop**: resolve this candidate's `(peptide_idx, charge)` ->
   `(start, end)` via the dense index (reusing `iim_dense_slot`), then
   unpack the sparse `(annotation_id, intensity)` pairs at that range into
   a stack-allocated `predicted_dense: Option<[f32; 174]>` (`None` if not
   configured, or no cache entry for this candidate — coverage gaps are
   silent, not errors). `observed_dense: Option<[f32; 174]>` is only
   allocated at all when `predicted_dense` is `Some` — no point tracking
   observed intensities for a candidate with nothing to compare against.
2. **During the loop**: for every theoretical fragment SAGE already tries
   to match against an observed peak (the pre-existing
   `select_most_intense_peak` binary-search call — **no new peak lookup is
   added**, this reuses the same match the existing hyperscore accounting
   already found), if a peak matched, write `observed_dense[annotation_id]
   = peak.intensity` at the same slot the prediction lives at. Predicted
   and observed are never separately "matched" against each other — they
   land in the same array index because both are derived from the same
   `(kind, idx, charge)` triple in the same loop iteration.
3. **After the loop**: `ms2_similarity::entropy_similarity(&observed,
   &predicted)` once, stored on `Score` (a plain field, not part of
   `combined_score`/ranking — `Score::rank_key`/`Ord` are untouched).

Both dense arrays are fixed-size stack arrays (696 bytes each), not heap
allocations — created and dropped within one `score_candidate` call, no
accumulation across candidates/spectra/threads. The one real shared
resource is `arrays.mmappet` itself, mapped exactly once by `Runner` and
shared read-only via `&'static` slices, same sharing model as
`IndexedDatabase` and the RT/IIM dense arrays.

## `ms2_similarity::entropy_similarity` (`crates/sage/src/ms2_similarity.rs`)

Unweighted spectral entropy similarity (Li et al. 2021, *Spectral entropy
outperforms MS/MS dot product similarity for small-molecule compound
identification*, Nature Methods) — `1 - JSD(P, Q) / ln(2)`, where `P`/`Q`
are `observed`/`predicted` each normalized to sum to 1 and `JSD` is their
Jensen-Shannon divergence (`H((P+Q)/2) - 0.5*H(P) - 0.5*H(Q)`). Range
`[0, 1]`; `0.0` when either vector sums to zero (an empty comparison has no
meaningful similarity, not perfect similarity) rather than `NaN`/`1.0`. A
slot where both vectors are `0.0` (fragment position beyond this peptide's
real length, or outside the model's vocabulary) contributes nothing to any
of the three entropies — safe to leave both dense arrays zero-initialized
for slots a given peptide never visits, no separate "valid slot" mask
needed. Unit-tested directly: identical vectors -> `1.0`, disjoint support
-> `0.0`, scale invariance (post-normalization), all-zero vector -> `0.0`
not `1.0`, partial overlap strictly between the two extremes.

## `Feature::ms2_entropy_similarity`: output wiring

New field, `0.0` default (same convention as `predicted_rt_external`/
`delta_rt_z2_external`), appended in three places kept in sync manually
(no shared macro across these three formats in this codebase):

- `results.sage.tsv` (`serialize_feature`/`write_features`,
  `sage-cli/src/runner.rs`) — appended at the very end of the column list.
- `results.sage.pin` (`serialize_pin`/`write_pin`, same file) — inserted
  right before `Peptide`/`Proteins`, which must stay the last two columns
  (Percolator/mokapot convention) — **not** appended at the end like the
  TSV.
- parquet (`sage-cloudpath/src/parquet.rs`'s `build_schema`/`write_col!`
  macro) — added right before the optional `reporter_ion_intensity` list
  group, matching `Feature` struct field order.

**Not** added to the LDA (`ml/linear_discriminant.rs`) or to
`combined_score`/`RankingScore` — this field is reported only, per the
explicit "no hard filtering on those" decision this feature was scoped
under. If it's ever worth feeding into ranking or the LDA, that's a
separate, deliberate follow-up, not implied by this change.

## Verified end-to-end, real F9477 data (2026-08-31)

Full workspace test suite green (`cargo test --workspace --all-features`) —
2 new `sage-cloudpath` tests (`fragment_intensity_cache`), 7 new `sage-core`
tests (`ms2_similarity`, including the annotation-id kind-block boundary
and every entropy-similarity edge case above), all 14 pre-existing
`sage-core`/`sage-cli` integration tests unaffected (their `Scorer` literals
just needed the three new `None` fields added).

Real search, full F9477 (same `pmsms`/`precursors`/`sage_config` as an
existing, previously-verified `run_sage` node — comparable baseline: 141.3s
without this feature configured at all):

- **206,991** PSMs, **206,197** (99.6%) with a nonzero
  `ms2_entropy_similarity` (the rest: peptides too long for the model, or
  fragment charges/kinds outside its vocabulary).
- Range `[0.0, 0.962]`, mean `0.374`, median `0.331` — sane, no `NaN`/`inf`,
  not degenerate (not all-0 or all-1).
- **Real discriminating signal, not noise**: mean `0.421` for targets
  (`label=1`, n=121,999) vs. `0.305` for decoys (`label=-1`, n=84,992) —
  genuine target/decoy separation in the expected direction, a first sanity
  check that this is measuring something real before any downstream
  consumer (mokapot/`sagepy-rescore`) is built on it.
- Runtime: 159s (TSV run) / 157s (PIN run) vs. 141.3s baseline — **+17.7s,
  ~12.5%**. Not yet profiled *where* that goes (loading/resolving the
  16.87M-entry index vs. per-candidate work) — the per-candidate cost
  itself should be O(1) dense-array reads plus a small fixed-size unpack,
  the same cost class as existing IIM eviction, so the one-time
  index-load/HashMap-construction step is the more likely dominant cost;
  not confirmed by profiling, flagged as a follow-up if this overhead ever
  matters in practice.
- `results.sage.pin`'s column order double-checked directly: `...,
  posterior_error, ms2_entropy_similarity, Peptide, Proteins` — confirmed
  the required-last-two-columns convention is preserved.

## Full MSBooster parity (2026-08-31, second pass)

After the first slice (`ms2_entropy_similarity` only) shipped, the user
asked for the rest of MSBooster's MS2-similarity feature set. Two research
passes against the real source (`github.com/Nesvilab/MSBooster`,
`src/main/java/features/spectra/SpectrumComparison.java`) were needed — the
first (pre-compaction, summarized rather than re-verified) wrongly claimed
MSBooster has *no* cosine/dot-product family at all; a second, source-level
pass (`gh api`/`raw.githubusercontent.com`, full 1011-line file read)
corrected this and pulled exact formulas, corcting the record.

**12 new `Feature` fields**, all computed from the same compacted
real-fragment-positions-only `(observed_real, predicted_real)` pair
`ms2_entropy_similarity` already used (`ms2_top6_matched_intensity`
additionally needs the full raw spectrum, see below) — no new `Scorer`
config, no new CLI flags. (Originally these took the full `[f32; 174]`
dense arrays directly; changed after a real bug was found -- see "Bug
found and fixed" below.)

| Field | MSBooster function | Notes |
|---|---|---|
| `ms2_weighted_entropy_similarity` | `weightedSpectralEntropy` | `entropy_similarity * H(one_normalize(predicted))^0.5` — **not** m/z-frequency-weighted (see below), a self-weighting transform. Unbounded above 1.0 by design. |
| `ms2_heuristic_entropy_similarity` | `heuristicSpectralEntropy` | If predicted-spectrum entropy `< 1.75`, reweight both vectors by `intensity^(H/2.75)` before recomputing entropy; else identical to `ms2_entropy_similarity`. |
| `ms2_cosine_similarity` | `cosineSimilarity` | Raw (non-normalized) vectors. |
| `ms2_dot_product` | `dotProduct` | Cosine of L2-normalized vectors — algebraically equal to `ms2_cosine_similarity` when both are defined; verified identical in the real run below. |
| `ms2_spectral_contrast_angle` | `spectralContrastAngle` | `1 - (2/π)·acos(cosine)`. |
| `ms2_euclidean_similarity` | `euclideanDistance` | Despite the MSBooster name, already a similarity (`1 - ||Δ||₂` on unit vectors) — can go negative (min `1-√2`). |
| `ms2_bray_curtis_similarity` | `brayCurtis` | `1 - Σ|Δ|/Σsum`, on unit vectors. |
| `ms2_pearson_corr` | `pearsonCorr` | **Real-fragment-positions only**, not the full 174-slot array (see below) — `-1.0` sentinel if degenerate. |
| `ms2_spearman_corr` | `spearmanCorr` | Same restriction; average-rank tie handling. |
| `ms2_hypergeometric_probability` | `hypergeometricProbability` | Scope-restricted (see below); `-log10(P(X >= k))`, no external stats crate needed (exact log-space PMF summation, `ms2_similarity::hypergeometric_upper_tail`). |
| `ms2_intersection` | `intersection` | Scope-restricted; raw count (not Jaccard), `u32`. |
| `ms2_top6_matched_intensity` | `top6matchedIntensity` | Numerator scope-restricted (top-6 *real* positions by predicted intensity); denominator uses the **full raw spectrum** (`query.peaks`, not just matched-fragment positions) -- the one metric here that needs more than the two 174-slot dense arrays. |

### Three deliberate deviations from MSBooster's literal definitions

1. **Pearson/Spearman restricted to real positions.** MSBooster's own
   vector is exactly as long as the model's real prediction count for that
   peptide — no padding. This repo's dense arrays are fixed at 174 slots
   with zero-padding for positions beyond a peptide's real length. Every
   *sum-based* metric above (cosine, dot product, Euclidean, Bray-Curtis,
   entropy) is mathematically unaffected by that padding — a `(0, 0)` slot
   contributes exactly `0` to every sum term. Pearson/Spearman are **not**
   sum-invariant to padding (they normalize by sample size via mean/
   variance/rank) — extra `(0, 0)` "agreement" points would silently bias
   the correlation upward. Fixed by tracking `is_real_dense: [bool; 174]`
   in `score_candidate` (set unconditionally per `(idx, charge)` the
   fragment loop visits, regardless of match) and filtering both vectors to
   just those positions before calling `pearson_corr`/`spearman_corr`. Real
   run: **zero** `-1.0` sentinel triggers across 206,991 PSMs — real
   peptides always had enough real positions.

2. **Hypergeometric probability/intersection use a narrower "population".**
   MSBooster's population is every theoretical fragment across all 6 ion
   kinds within the scan's observable m/z range — a separate, larger
   computation than what the job's own search uses. Per explicit decision
   (avoid generating fragment kinds beyond what's actually configured),
   both are restricted to this repo's 174-slot vocabulary + `is_real`:
   population = real positions, population successes = how many were
   observed, sample = the subset the model actually predicted for (nonzero
   `predicted`), sample successes = how many of *those* were observed. This
   is a real, if narrower, coherent hypergeometric test — not degenerate,
   since Prosit `compact_trt`'s own cache coverage is itself genuinely
   sparse (see `fragment_slots.py`'s `ordinal_in_range` finding: roughly
   half of a peptide's structurally-valid slots may still lack a cache
   entry), so sample size is meaningfully smaller than population in
   practice.

3. **`getWeights`/`freqs`-dependent `weighted*` variants are skipped
   entirely** (`weightedCosineSimilarity`, `weightedDotProduct`,
   `weightedSpectralContrastAngle`, `weightedEuclideanDistance`,
   `weightedBrayCurtis`, `weightedPearsonCorr`) — these multiply in an
   m/z-binned frequency table whose own computation wasn't found in either
   research pass (not in `SpectrumComparison.java`). Shipping a guessed
   weighting scheme would be silently wrong, not just approximate, so
   these are left out rather than faked. `weightedSpectralEntropy`/
   `heuristicSpectralEntropy` are **not** part of this skip — they
   self-weight from their own computed entropy, no `freqs` dependency.

### Verified end-to-end, real F9477 data (2026-08-31)

Full workspace test suite green: 36 `ms2_similarity` unit tests (up from 7
— one genuine test bug caught and fixed during this pass: `intersection`
with `top_n >= population` correctly includes *every* real position
regardless of its own observed value, which the first draft of
`intersection_zero_without_overlap` didn't account for).

Real search, same F9477 `pmsms`/`precursors`/`sage_config` as the
`ms2_entropy_similarity`-only run: **159s for the first 11 new metrics —
identical to the single-metric run**, no measurable additional overhead
(all cheap scalar arithmetic over already-resident `[f32; 174]` arrays,
plus a small bounded `Vec` for the Pearson/Spearman real-position subset).
Adding `ms2_top6_matched_intensity` afterward (needs the full raw
`query.peaks` list, not just the two dense arrays) brought it to **163s**
— a small, real ~2.5% bump from the per-candidate peak-intensity
collection + sort, still modest.

**Every metric shows the same real target/decoy separation direction**
(206,991 PSMs, label=1 vs label=-1 means): cosine `0.384` vs `0.249`,
spectral contrast angle `0.273` vs `0.166`, Bray-Curtis `0.309` vs
`0.203`, Pearson `0.249` vs `0.078`, Spearman `0.272` vs `0.129`,
hypergeometric `0.692` vs `0.430`, Euclidean `-0.063` vs `-0.210` (still
target > decoy despite both being negative), intersection `13.84` vs
`13.65` (real but much smaller separation, expected — a top-20 count is
naturally close to saturated for both). `ms2_cosine_similarity` and
`ms2_dot_product` matched to float precision in every sampled row, as
expected algebraically. `ms2_heuristic_entropy_similarity` matched
`ms2_entropy_similarity` exactly for every sampled high-hyperscore
candidate (predicted-spectrum entropy `>= 1.75` for these, so the
reweight branch never triggers) — plausible, not yet checked against a
real low-entropy (short/simple predicted spectrum) case specifically.
`ms2_weighted_entropy_similarity` exceeded `1.0` in real data (max
`1.566`), confirming the intentionally-unbounded behavior.

## Verified together with `--predicted-rt` on the real full second pass (2026-08-31)

All runs above used the 250,000-row `select_recalibration_precursors` subset
(the confident-hit-selection stage of this project's normal 2-pass
pipeline) — the same input the pre-existing baseline node used, so the
timing comparisons were fair, but not the largest real search this project
runs. Re-verified against the actual **second/final pass**: the real
`run_sage` node (`cf8172b0f7c258bf...`) that already has `--predicted-rt`
wired in (RT hard eviction via `evict_rt_iim_mismatches` + `delta_rt_z2_external`
in the LDA/`combined_score`, see `predicted_rt_iim.md`) — full
`rt_corrected_precursors.mmappet` (not a confident-hit subset),
`recalibrate_pmsms_mz`'s pmsms, `--annotate-matches --write-pin`, real
`predicted_rt.parquet`, plus this feature's two new flags added on top of
that exact real command (pulled verbatim from the node's own
`.rip/dependencies.toml`).

- Baseline (this exact node, before this feature existed): **376.8s**
  (necroflow's own recorded `run.toml`).
- With `--predicted-fragment-intensity-*` added: **412s** (+35.2s, ~9.3%).
- **468,368 PSMs** — the real full precursor set, not the 250K-row subset.
- `predicted_rt_external` populated for 100% of rows (RT feature/eviction
  active throughout), `ms2_entropy_similarity` for 99.65%.
- `delta_rt_z2_external` correctly lower for targets (`1.209`) than decoys
  (`1.407`) — RT hard eviction behaving as already established.
- **Every one of the 12 MS2 metrics again shows the correct target/decoy
  separation direction, and noticeably *stronger* than the
  250K-precursor-subset run** (e.g. entropy similarity `0.478` vs `0.319`
  here, vs. `0.421` vs `0.305` there) — consistent with RT hard eviction
  having already filtered out weaker candidates before MS2 similarity is
  even computed. Confirms the two features compose correctly and don't
  interfere with each other (RT eviction still runs on `InitialHits` before
  `score_candidate`'s MS2 unpacking; the MS2 metrics never see candidates
  RT already rejected).

## Bug found and fixed: predicted covers charges the job never searches (2026-08-31)

Caught by a direct question, not by testing: **is `predicted`/`observed`
compared against `is_real` before summing, for every metric, or only
some?** The answer was only some -- and the ones that weren't had a real,
currently-active bug.

The cache's `predicted_dense` is unpacked directly from
`arrays.mmappet`'s sparse entries, which cover all 3 Prosit fragment
charges (1, 2, 3) for a given `(sequence, precursor_charge)` key,
independent of what *this job's own* `max_fragment_charge` config
actually searches. The real F9477 production config uses
`max_fragment_charge: 1` -- `fn max_fragment_charge` (`scoring.rs`) maps
that to a loop bound of `1..2`, so SAGE's own fragment-generation loop
only ever visits **fragment charge 1** (roughly 1/3 of the 174-slot
vocabulary) for every real job run so far. `is_real_dense` (correctly)
only marks those charge-1 slots true. But `predicted_dense` -- unpacked
straight from the cache, with no awareness of `max_fragment_charge` at
all -- has real, nonzero Prosit predictions at the charge-2/3 slots too.

Eight of the twelve MSBooster-parity metrics (`entropy_similarity`,
`weighted_entropy_similarity`, `heuristic_entropy_similarity`,
`cosine_similarity`, `dot_product`, `spectral_contrast_angle`,
`euclidean_similarity`, `bray_curtis_similarity`) summed over the **full**
174-slot dense array with no `is_real` filtering at all -- relying on an
invariant ("padding with `(0, 0)` pairs is neutral for every sum-based
metric") that's only true when `predicted` is actually `0` at every
non-real slot. It wasn't: those charge-2/3 slots had real `predicted > 0`
values silently compared against phantom `observed = 0` (SAGE never even
attempts to match a fragment charge it doesn't search), biasing all eight
metrics downward for every real job run to date. The other five
(`hypergeometric_probability`, `intersection`, `top6_matched_intensity`,
`pearson_corr`, `spearman_corr`) were **never** affected -- each already
took an explicit `is_real` parameter and filtered through it before ever
reading `predicted`.

**Fix**: rather than adding an `is_real` parameter to each of the eight
(a "remember to mask" pattern that already failed once), every metric
function in `ms2_similarity.rs` was changed to take plain `&[f32]` slices
instead of `&[f32; N_FRAGMENT_SLOTS]` arrays, and `score_candidate` now
builds **one** compacted `(observed_real, predicted_real)` pair (the same
compaction Pearson/Spearman already did) and feeds it to all thirteen
metrics uniformly -- there's no code path left where a metric can see an
unfiltered slot at all, so this class of bug can't recur by omission.
`hypergeometric_probability`/`intersection`/`top6_matched_intensity`
dropped their now-redundant `is_real` parameter (population membership is
simply `observed_real.len()` once the slice is already compacted).

Real measured effect (same F9477 second-pass run as above, target/decoy
means): `entropy_similarity` `0.478/0.319` → `0.533/0.363`,
`cosine_similarity` `0.453/0.267` → `0.484/0.291`, `euclidean_similarity`
`0.008/-0.194` → `0.044/-0.172`, `bray_curtis_similarity`
`0.361/0.215` → `0.404/0.247` -- all shifted up ~10-15%, as expected once
phantom mismatches are excluded. The other five metrics are
byte-for-byte identical before and after, confirming they were never
exposed. Runtime: 423s (vs. 412s before this fix -- +2.7%, within
run-to-run noise for this machine).

## Explicitly out of scope here

- The `getWeights`/`freqs`-dependent `weighted*` variants (see above) —
  would need a third research pass to find where MSBooster computes its
  m/z-binned frequency table before these can be faithfully ported.
- Any downstream consumer (mokapot rescoring, `sagepy-rescore`) actually
  using any of these 13 fields — same "reported but unconsumed so far"
  status `predicted_rt_external`/`delta_rt_z2_external` had when they
  first landed.
- Profiling/optimizing runtime overhead — none measured so far beyond
  wall-clock comparison of three runs (141.3s baseline, 159s with the
  first 11 new metrics, 163s with `ms2_top6_matched_intensity` added).
