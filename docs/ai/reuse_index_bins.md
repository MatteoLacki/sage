# Reusing fragment-index pages across queries (`experimental/reuse_index_bins`)

Two related optimizations to `IndexedDatabase::query`/`page_search`'s
per-query cost, on the two-level bucket index described in
`docs/ai/dump_peptides.md`'s sibling discussion (see also
`docs/ai/predicted_fragment_intensity.md`). Stage 1 is implemented and
verified on this branch. Stage 2 is a design, not yet built.

## Background: what's redundant today, even after Stage 1

`page_search(mass, tol)` (`database.rs`) does two binary searches per
call: an **outer** one over `min_value` (page boundaries — depends only on
`(fragment_mass, tol)`) and an **inner** one within the selected page
(depends only on `(page, pre_idx_lo, pre_idx_hi)` — the spectrum's fixed
precursor window). Called once per `(peak, fragment_charge)` in
`matched_peaks_with_isotope` (`scoring.rs`), both searches were redone
from scratch on every call, even when consecutive peaks (spectrum peaks
are pre-sorted by mass, `spectrum.rs`) land on the same page — likely,
since bucket width at this pipeline's `bucket_size` often exceeds the
fragment tolerance width (see the `bucket_size`/L1-L2 discussion this
design grew out of).

## Stage 1 (implemented, this branch): reuse within one spectrum

`IndexedQuery::page_search_batch(&self, windows: &[(f32, Tolerance)])`
takes every `(mass, fragment_tol)` window for one spectrum (all peaks x
all fragment charges) at once instead of one call per window:

1. Resolve each window's page range via the outer search (unchanged
   per-window cost here — see Stage 2 for removing this too).
2. Collect `(page, window_index)` pairs, sort by page.
3. Walk pages in order; for each *unique* page, do the inner
   (precursor-scoped) binary search **once**, then apply every window
   referencing that page against the same resolved `[inner_left,
   inner_right)` slice.

`pre_idx_lo`/`pre_idx_hi` are fixed for the whole `IndexedQuery` (one per
`(spectrum, isotope_error)` call to `matched_peaks_with_isotope`), so this
is exact, not approximate — same results, less redundant work.

**Verified, and the headline number was misleading (2026-09-08).** All 131
`sage-core` unit tests + 14 integration tests pass unchanged, and output
is byte-identical to the pre-Stage-1 baseline at full scale (363,156 PSMs,
identical id set, all 54 columns byte-identical) — correctness is solid.
Performance is not a win: the first benchmark (10K precursors, "top 10K by
intensity" subset — matching `quick_test_10k.toml`'s own selection) showed
search phase 9.30s → 4.73s (~49% faster), but that subset is biased.
Rerun on a **random** 10K-precursor sample and on the **full**
716,614-precursor set (both cleanly isolated — identical binary lineage,
identical input files, only the Stage 1 code differing):

| subset | search, no Stage 1 | search, +Stage 1 | delta |
|---|---|---|---|
| top-10K-by-intensity | 9.30s | 4.73s | −49% |
| random 10K | 4.31s | 4.49s | **+4.1%** |
| full 716,614 | 304.76s | 317.51s | **+4.2%** |

Random and full-scale agree with each other and disagree with the biased
sample. Root cause (plausible, not separately profiled): top-intensity
spectra have systematically more peaks (cleaner signal passes
`min_peaks`/`max_peaks` more readily), so they're exactly the case where
multiple peaks land on the same page — real sharing opportunity. Typical
spectra have fewer peaks, little to share, and `page_search_batch`'s fixed
per-call cost (three heap allocations — `bounds`, `page_window`,
`results` — plus a sort, versus the original `page_search`'s
zero-allocation lazy iterator chain) isn't recouped there.

**Follow-up fix (2026-09-08, same day): stop allocating fresh `Vec`s per
call.** The first implementation built three owned `Vec`s per
`(spectrum, isotope_error)` call (`bounds`, `page_window`, and a
`results: Vec<&Theoretical>` the caller then iterated) — a first
back-of-envelope check (up to ~2.8M calls in a full run, each a few
small allocations) suggested pure allocation cost was too small to fully
explain the regression on its own, pointing at either allocator
contention under 16 concurrent threads, or just the sort/bookkeeping cost
not being recouped when a spectrum's peaks don't actually share pages
(the common case). Reworked to: thread-local scratch buffers for `bounds`
and `page_window` (rayon's worker threads are long-lived OS threads, so
`thread_local!` persists correctly across the many spectra one worker
processes — no allocation on repeat calls, just `.clear()` and reuse), and
a callback (`on_match: impl FnMut(&Theoretical)`) instead of returning a
`Vec` at all, restoring the same zero-allocation-for-results shape
`page_search`'s original lazy iterator chain had. `scoring.rs`'s own
per-spectrum `windows: Vec<(f32, Tolerance)>` got the same treatment.

| subset | positional-only | v1 (fresh `Vec`s) | v2 (scratch + callback) |
|---|---|---|---|
| top-10K-by-intensity | 9.30s | 4.73s (−49%) | 4.70s (−49%) |
| random 10K | 4.31s | 4.49s (+4.1%) | 4.38s (**+1.6%**) |
| full 716,614 | 304.76s | 317.51s (+4.2%) | 306.68s (**+0.63%**) |

Removing the allocations closed nearly the entire regression (full-scale:
+4.2% → +0.63%, essentially noise) while leaving the genuine win on
many-peak spectra untouched, as expected — that win was always about
avoiding redundant inner binary searches, never about allocation. Full
716,614-precursor correctness reverified after this change too: 363,156
PSMs, identical id set, all 54 columns byte-identical to the
pre-Stage-1 baseline. All 131 unit + 14 integration tests pass.

**Conclusion, revised: v2 looks worth merging.** Near break-even
(≤1.6%, likely partly the residual `sort_unstable_by_key` cost, not
separately profiled) on both realistic samples tested, a real ~49% win
when spectra have enough peaks for genuine page-sharing, correctness
fully verified, and the allocation-removal lesson is worth keeping in
mind for Stage 2's design (batch across spectra) from the start, rather
than re-discovering it there too.

## Stage 2 (design, not implemented): merge across spectra

Stage 1 only shares work *within* one spectrum. The outer search (page
resolution by fragment mass) doesn't depend on precursor at all — it's
shareable across *any* spectra whose peaks happen to land near each other
in fragment-mass space, which across a whole run's worth of spectra is
common in dense mass regions. The inner search depends on precursor mass
too, but that's shareable across spectra with *similar* precursor mass —
which is exactly what batching by precursor-mass locality can arrange for.

### Data model

Extend the query unit from "one spectrum" to "one batch of `B` spectra",
where a spectrum entry is really `(spectrum_id, isotope_error)` — since
each `isotope_error` value shifts the effective precursor mass and thus
`pre_idx_lo`/`pre_idx_hi`, it needs its own precursor window, but **not**
its own peak list (fragment masses don't depend on isotope error — reuse
one merged peak stream per spectrum across all its isotope-error variants,
only re-deriving `pre_idx_lo`/`pre_idx_hi` per variant). `chimera`/
`wide_window` search modes use a different code path
(`score_chimera_fast`, `scoring.rs`) — out of scope here, not silently
broken, just not addressed by this design; revisit separately if wanted.

### Batching: sort by precursor mass, not file order

Sort every `(spectrum, isotope_error)` entry in the run by its *effective*
precursor mass (`precursor_mass - isotope_error * NEUTRON`) once, up
front — cheap (`O(n log n)` over spectrum count, negligible next to the
search itself). Chunk that sorted sequence into batches of size `B`. This
is the key move for sharing the *inner* search: spectra within one batch
now have precursor windows that are close together (often overlapping),
not arbitrary.

`B` is a real tuning knob: bigger batches mean more sharing opportunity
(more spectra's windows landing on the same page) but coarser-grained,
more serialized per-batch work and higher per-batch latency. Pick `B` so
`num_batches` stays well above `num_threads` (e.g. `B` such that
`num_batches >= 4-8x` core count) for load balancing — start there and
tune empirically against real F9477 data, same as `bucket_size` was.

### The merge itself

For one batch:

1. **k-way merge the batch's peak streams** into one globally-sorted
   `(fragment_mass, tol, batch_local_id)` stream. Each spectrum's own
   peaks are already sorted, so this is a textbook k-way merge (`k = B`)
   — a binary heap of per-spectrum iterators for larger `B`, or plain
   collect-and-sort for small `B` where the constant-factor difference
   won't matter. Isotope-error variants of the same spectrum reuse the
   same peak entries, tagged additionally with which variant(s) apply.
2. **One forward sweep against `min_value`** resolves every window's page
   range via a monotonic two-pointer walk (interval join, not per-window
   binary search) — `O(windows + pages)` for the *whole batch*, replacing
   Stage 1's still-per-window outer search.
3. **Group by page**, same as Stage 1. Within a page, further group by
   *distinct* `(pre_idx_lo, pre_idx_hi)` among the batch's entries
   touching that page (adjacent/overlapping windows from precursor-mass
   locality should often collapse to few distinct ranges) — one inner
   binary search per distinct range instead of per `(spectrum,
   isotope_error)` pair.
4. **Scatter final matches** back to each `(spectrum, isotope_error)`'s own
   `InitialHits`/`PreScore` accumulator via the tag carried through the
   merge — same bookkeeping shape as today, just fed by a shared
   resolution pass instead of `B` independent ones.

### Parallelism: the real complication

Today's parallel axis is "one thread = one spectrum"
(`search_processed_spectra`, `runner.rs`). Stage 2's merge couples a whole
batch's work together, so the axis has to move to "one thread = one
batch" — `par_iter()` over `chunks(spectra_sorted_by_precursor_mass, B)`
instead of over individual spectra. This is a natural, consistent
generalization (today's per-spectrum work is already single-threaded
within its own thread; this just grows the single-threaded unit from 1
spectrum to `B`), but it's the part most likely to hide a bad surprise —
a `B` that's too large could starve threads late in a run (fewer, coarser
batches finishing at uneven times), and the sort-by-precursor-mass
preprocessing step touches every spectrum in the run, so its own cost and
memory footprint need checking at real scale before trusting the design.

### Build order

1. Land Stage 1 on its own (already done) — confirm full-scale
   correctness (in progress) and keep it as an independent, revertible
   commit.
2. Prototype the k-way merge + page/inner-range sharing in isolation
   (unit tests with small synthetic multi-spectrum fixtures, comparing
   against today's per-spectrum-independent results) before touching the
   parallelism model at all — get the merge algorithm provably correct
   first.
3. Only then move the parallel axis from spectrum to batch, profiling `B`
   against real F9477 data (10K-precursor job first, matching this whole
   session's methodology, then full-scale byte-identical PSM
   verification) before trusting timing claims.
4. Explicitly decide what happens to `chimera`/`wide_window` search
   (separate code path, not addressed by this design) rather than let it
   silently diverge.

### Open questions to resolve during prototyping, not before

- Real sweet spot for `B` — no basis to guess it without real profiling.
- Whether k-way-merging via a binary heap is worth its complexity over
  plain collect-and-sort at realistic `B` — depends on how large `B` ends
  up being.
- Whether the precursor-mass sort should be global (whole run) or
  per-file/per-chunk — global gives more sharing opportunity but couples
  more of the pipeline together; per-file is more local/simpler to reason
  about and may capture most of the benefit if precursor masses are
  already reasonably spread within a file.
