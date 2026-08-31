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

## Explicitly out of scope here

- Any additional MSBooster/DiaTracer-style features (spectral-angle,
  intersection count, hypergeometric probability, top6matchedIntensity) —
  `ms2_entropy_similarity` is the first slice, not the whole set. The
  `predicted_dense`/`observed_dense` pair computed in `score_candidate` is
  general-purpose; adding another feature from the same two arrays is a
  small, localized change, not a redesign.
- Any downstream consumer (mokapot rescoring, `sagepy-rescore`) actually
  using this field — same "reported but unconsumed so far" status
  `predicted_rt_external`/`delta_rt_z2_external` had when they first
  landed.
- Profiling/optimizing the ~12.5% overhead above.
