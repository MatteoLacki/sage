# Fragment grid: precursor blocks × fragment cells (2026-10-16)

Necromerge2 plans: `plans/fragment_grid.md`, `plans/grid_match_simd.md` (the
branch-free match kernel), `plans/grid_compact_cells.md` (equal-count cells and
16-bit fragments). A second layout of the fragment index
(`crates/sage/src/fragment_grid.rs`), chosen when the index is built. Prototype:
switched by environment variables, default unchanged.

## Why

Sage's own layout (`database.rs`) is fragment-major: fragments sorted by m/z, cut
into `bucket_size` pages, each page sorted by peptide id, which is precursor mass
order. Every fragment window binary-searches the peptide-id column of each page it
touches for the precursor range, a chain of dependent DRAM probes
(`resolve_page_peptide_bounds`, `simd.md`). That layout suits wide precursor windows.
Ours are narrow (±6.7 ppm, isotope errors 0–3), so the grid turns it round:
precursor-major blocks, fragment cells inside, found by arithmetic.

## Layout

```
                 cells: runs of fine keys (mz.to_bits() >> (23 - f)), ~T fragments per block-cell
                 c0     c1     c2    ...
block 0  rows 0..R    [run]  [run]  [run]  ...
block 1  rows R..2R   [run]  [run]  [run]  ...
cell_start[block][c] -> fragment offset (relative to the block) into its 16-fragment groups
```

- **Blocks:** R consecutive peptide rows (R = 4096). Rows are mass-sorted, so a block
  is a precursor-mass band; equal row counts keep blocks equally full where peptide
  density varies. A precursor window reaches the grid only through rows: sage's binary
  search over the per-row masses gives the candidate rows, `row / R` the blocks.
- **Fine keys:** for non-negative floats the bit pattern is monotone in the value, and
  its top bits are the exponent and leading mantissa bits. `mz.to_bits() >> (23 - f)`
  keeps f mantissa bits: 2^f fine keys per octave, about constant ppm wide (f = 14:
  61–122 ppm).
- **Cells (equal-count):** consecutive fine keys grouped so each cell holds about T
  fragments per block on average over the grid (T = 8), from a histogram of all
  fragments' fine keys. A cell closes once it holds T × blocks fragments or spans
  2^(f−7) fine keys, so no cell is wider than 2^16 bit patterns. `cell_of_fine` maps
  each fine key from the lightest fragment's to the heaviest's to its cell (`u16`,
  ~115 KB for class I); `cell_first_fine` holds each cell's first key. Window keys are
  clamped to that range, and a window whose high end is not positive holds nothing.
  Cells are contiguous fine-key ranges, so a window's cells are
  `cell_of_fine[fine(lo)] ..= cell_of_fine[fine(hi)]`, and any such partition is exact;
  the counts only decide how much is scanned.
- **Inside a cell:** generation order, i.e. peptide order. Nothing relies on it.
- **Storage:** `FragmentGrid` owns the fragments as `FragmentGroup`s, 64 B each (one
  cache line, `align(64)`): 16 m/z residuals (`u16`), then their 16 rows within the
  block (`u16`, R ≤ 65,536). A residual is `mz.to_bits() − cell_first_bits(cell)`
  (below 2^16 by the width cap), so `cell_first_bits + residual` is the exact f32, and
  within a cell residual order is m/z order. Each block starts on a group; lanes past
  its last fragment are never read. Also block group starts and per-block `cells + 1`
  cell starts (`u32` fragment offsets relative to the block). 4 B per fragment.
  `IndexedDatabase::grid` holds it; `fragment_mzs`/`fragment_peptide_indices` and
  `min_value` stay empty. `size()` and `fragments()` read the grid and rebuild the
  exact m/z and row (`database_size` and `dump_fragment_index` read the page columns,
  so they see nothing for a grid).

## Build

`Parameters::build_from_peptides` with `fragment_layout = Grid` calls `build_grid`.
The fragments come from `Parameters::indexed_fragments`, shared with the page layout,
so their m/z bits are identical. `FragmentGrid::build` never holds more than the grid:
1. **Count:** per block, in parallel, generate its fragments into thread-local
   scratch; keep the count and the lowest/highest fine key, and add every fine key to a
   histogram (rayon fold/reduce of per-split histograms over a fixed 16–32,768 m/z
   range; a fragment outside it is a hard error).
2. **Cells and offsets:** `equal_count_cells` over the histogram; block group starts
   by prefix sum.
3. **Place:** the groups are allocated (a parallel collect); per block, in parallel,
   into its own disjoint groups: generate again, count per cell, prefix sum (the
   block's cell starts), place in generation order (a stable counting sort): residual
   and local row into lane `offset % 16` of group `offset / 16`.

No global sort. The page layout instead sorts every fragment by m/z, then sorts each
page by peptide, and peaks at about 16 B per fragment for those copies.

## Search

`matched_peaks_with_isotope` asks `IndexedQuery::matched_candidates` for the matched
candidates' offsets from `candidate_rows().start` (the exact rows
`[exact_pre_idx_lo, pre_idx_hi)`), one entry per (window, fragment) match, and counts
them into `preliminary` (`InitialHits::add_matches`, below). The grid answers with
`FragmentGrid::matched_candidates`; pages run their own `page_search_batch` with a
closure that pushes offsets. Per call:
- per window, once (`cell_windows`): `(lo, hi)`, its cells `first..=last`, and its
  bounds as residuals: `residual_lo` in the first cell (0 when `lo` lies below it),
  `residual_hi` in the last (`u16::MAX` when `hi` lies above it);
- per block overlapping the candidate rows (outer loop): the candidate rows local to
  the block, `[a, b)`, and the offset that turns a local row into a candidate offset;
- per window (inner loop), `CellWindow::runs`: the first cell's fragments need
  `residual >= residual_lo`, the last cell's `residual <= residual_hi` (both when the
  cells coincide), and the cells between lie inside the window by construction;
- each run as whole groups, the run's first and last groups masked by position
  (`lanes_in`), so another cell's residuals (not comparable) and lanes past the block
  are never used.

**Exact by construction.** Every fragment in the window lies in its cells (monotone
keys). Within one cell, `bits(lo) <= bits(m/z) <= bits(hi)` is
`residual_lo <= residual <= residual_hi` (positive floats order as their bits). The
row check is the candidate-range check. So the kernel selects exactly the (window,
fragment) pairs the pages select; only their order differs, and counting doesn't depend
on it.

`matched_candidates_avx2` (`#[target_feature(enable = "avx2,popcnt")]`, chosen at run
time; `matched_candidates_scalar` otherwise) per 16 fragments, without branches:
- residual mask: `vpmaxuw`/`vpminuw` + `vpcmpeqw` for both bounds, AND;
- row mask: subtract `a`, `vpminuw`/`vpcmpeqw` against `b − a − 1`;
- AND, `vpacksswb` with zero then `vpmovmskb` (lanes 0–7 in bits 0–7, lanes 8–15 in
  bits 16–23), AND with the position mask;
- per half: rows widened (`vpmovzxwd`) plus the block's offset, `vpermd` through
  `LEFT_PACK[mask]` (256 × 8 lane indices, a `const fn`), unaligned store into the
  matched list, advance by `popcnt` (`reserve` per run leaves room for 16 per group).
  No `pext` (microcoded on Zen 2).

`FragmentGrid::search` (the per-fragment callback, reached through
`page_search_batch` and `page_search`) rebuilds each fragment of the window's cells
exactly and compares the float m/z; the tests use it as the reference.

## Branch-free match handling (2026-10-16, `c03cecb`)

This section records the first branch-free kernel, over 8-lane groups of f32 m/z and
u32 rows with equal-ppm cells; the 16-lane residual kernel above replaced it.

Profiled first with in-source cycle counters (`perf` needs root): the grid kernel of
`062a871` was ~67% of class I search CPU, and re-running each call warm three ways
split it into scanning 26% (AVX2 compare, but `scan_masses_avx2` not inlined, one
call per run, and a scalar 0–7 tail with unpredictable branches), per-match
iteration 42% (`tzcnt` loop, a load from the separate id column, an unpredictable
edge-block id check), counting with a first-touch branch 20%, and the callback's
own overhead 12%. 268 G fragments scanned, 21.7 G matches.

`matched_candidates_avx2` (`#[target_feature(enable = "avx2,popcnt")]`, chosen at run
time; `matched_candidates_scalar` otherwise) handles a group in 19 instructions, the
loop branch the only branch:
- m/z mask: aligned load, `vcmpps >= lo`, `vcmpps <= hi`, AND;
- candidate mask: aligned load of the rows, subtract the first candidate row, unsigned
  `vpminud`/`vpcmpeqd` against the candidate count (the difference is the offset);
- AND, `vmovmskps`, `popcnt`;
- left-pack: `LEFT_PACK[mask]` (256 × 8 lane indices, built by a `const fn`),
  `vpermd` of the offsets, unaligned store into the matched list, advance by the
  popcount (`reserve` per run leaves room for 8 per group).

`add_matches` then updates `preliminary[first_slot + offset]` per entry without a
first-match branch (`scored_candidates += (matched == 0)`, the whole `PreScore`
rewritten). A per-candidate `u16` count array plus a pass over all candidates was
planned instead; the list avoids that pass (about as many candidates per call as
matches).

Without `popcnt` enabled the compiler emulated it (5 multiply/shift instructions per
group); `target_feature` on the whole kernel is what lets the comparison, packing and
loop inline into one loop.

### Verification

- **Unit:** `matched_candidates_equal_the_exact_search_in_every_lane_width`
  (`fragment_grid.rs`): random fragments with repeats, 46 windows, 8 candidate ranges,
  R ∈ {1, 3, 8, 16, 64} × m ∈ {0, 4, 10}: dispatched, scalar and AVX2 against
  `search`. `assert_search_matches_brute_force` (`database.rs`) also checks
  `matched_candidates` for pages and every grid shape. `add_matches` against the
  per-match callback it replaced. Breaking the m/z comparison (`>=` to `>`), the
  candidate bound, the run end rounding, the scalar bound or the packed count each
  fails them.
- **Suites:** sage-core, sage-cli and sage-cloudpath, pages and grid (R 4096/m 10,
  R 3/m 4, R 1/m 0).
- **Data** (`temp/fragment_grid/simd/check_all.sh`): the small inputs (8 formats,
  monolithic and sliced, G23), F9477 runs 10 and 11 and class I, grid and pages:
  all identical to `062a871`'s pages output.

### Measurements

F9477 class I sliced, search seconds, `062a871` grid / this kernel:

| threads | dense slices 1+2 | light slice 3 | search | wall |
|---:|---:|---:|---:|---:|
| 16 | 82.3 / 79.3 (1.04x) | 18.0 / 16.4 (1.10x) | 100.3 / 95.7 (1.05x) | 147.8 / 146.3 |
| 8 | 110.1 / 81.5 (1.35x) | 26.6 / 24.5 (1.09x) | 136.7 / 106.0 (1.29x) | 192.0 / 158.0 |
| 4 | 190.6 / 122.2 (1.56x) | 46.7 / 38.6 (1.21x) | 237.3 / 160.8 (1.48x) | 318.7 / 237.1 |

The dense slices are now bound by memory traffic at 16 threads: this kernel gains
nothing from 8 to 16 threads there (81.5 → 79.3 s), while the old one still did
(110.1 → 82.3 s, SMT hiding its dependent loads). The groups read 8 B per scanned
fragment instead of 4, which offsets most of the instructions saved; with fewer
threads, where bandwidth is not the limit, the new kernel is 1.35–1.56x faster. The
next lever is fewer bytes scanned, not fewer instructions. Build: 19.1 s against
17.7 s at 16 threads (the pad-filled allocation pass).

F9477 production search (monolithic tryptic, 1,539 blocks): run 11 112.2 s → 94.2 s
(pages, with the new counting: 125.0 → 122.2 s); run 10 with pages 23.3 → 22.9 s warm.

## Query order (sliced search)

Inside a slice, the runs of one spectrum's charges are sorted by their first query's
`window_low` (`QueryOrder::Mass`), so neighbouring searches read the same part of the
index. Results are kept by query index, so output cannot change. Measured: 5% faster
with pages, no difference with the grid (below).
Gzipped text inputs keep file order: a read out of order resumes from the last gzip
access point and re-inflates up to 4 MB (`IndexedSpectra::prefers_forward_reads`).
Before this change slices were searched in file order (chosen for gzip, never
revisited for the other formats).

## Switches (prototype)

| Variable | Values | Default |
|---|---|---|
| `SAGE_FRAGMENT_LAYOUT` | `pages`, `grid` | `pages` |
| `SAGE_GRID_ROWS` | R, at most 65,536 | 4096 |
| `SAGE_GRID_FINE_BITS` | f, 7–16 | 14 |
| `SAGE_GRID_CELL_FRAGMENTS` | T, at least 1 | 8 |
| `SAGE_QUERY_ORDER` | `mass`, `file` (sliced search) | `mass` |

`SAGE_GRID_CELL_BITS` (equal-ppm cells, up to `c03cecb`) is now an error.
`Builder::make_parameters` reads the layout into `Parameters::fragment_layout`,
`#[serde(skip)]`, so configs and `results.json` are unchanged. Each grid build logs
its shape and directory size (`fragment grid: ...`).

## Verification

- **Unit:** `grid_matches_brute_force_on_block_and_cell_edges` (`database.rs`): the
  page test's peptides, fragments and windows plus octave-edge fragments and windows,
  for R ∈ {1, 2, 3, 4, 6, 7} × m ∈ {0, 4, 10, 23}; databases without fragments.
  Breaking the run end, the edge-block id check, the cell clamp or the non-positive
  lower bound each fails it.
- **Suites with the grid forced:** `SAGE_FRAGMENT_LAYOUT=grid cargo test -p sage-core`
  (also R 3/m 4 and R 1/m 0) and `cargo test -p sage-cli` pass.
  `check_all_ions_visited` asserts the page structure only for pages.
- **Small inputs** (`temp/fragment_grid/check_small.sh`, necromerge2 root): the
  streaming test data's short cut in all 8 formats, monolithic and sliced, and G23
  tryptic monolithic and sliced: grid + mass order identical to pages + file order.
- **F9477 production** (regression runs 10 and 11: monolithic tryptic, file order; grid
  of 1,539 blocks × 4,700 cells, 28.9 MB directory): identical, 1,524,398 and 4,729,361
  matched fragments.
- **F9477 class I sliced** (unspecific 8–12 aa, 923,634 queries, 3 slices): every run
  below identical to pages + file order.

## Measurements (Ryzen 7 4800H, 8 cores / 16 threads, 62 GB)

F9477 class I sliced, totals over the 3 slices (`temp/fragment_grid/run_classI.sh`).
Slice boundaries move slightly between runs because the budget is half of the live
`MemAvailable`; results do not.

| layout, order | R | m | build | search | wall | peak anon / RSS |
|---|---:|---:|---:|---:|---:|---:|
| pages, file (before) | | | 58.9 s | 185.7 s | 280.6 s | 38.6 / 44.9 GB |
| pages, mass | | | 53.4 s | 177.2 s | 263.8 s | 39.2 / 49.7 GB |
| **grid, mass** | 4096 | 10 | **17.7 s** | **100.3 s** | **147.8 s** | **26.4 / 39.2 GB** |
| grid, file | 4096 | 10 | 17.9 s | 99.9 s | 148.2 s | 26.3 / 39.1 GB |
| grid, mass | 1024 | 10 | 17.7 s | 112.1 s | 159.5 s | 27.0 / 39.8 GB |
| grid, mass | 16384 | 10 | 19.2 s | 130.5 s | 180.1 s | 26.2 / 39.0 GB |
| grid, mass | 4096 | 8 | 18.0 s | 111.8 s | 158.4 s | 26.1 / 38.9 GB |
| grid, mass | 4096 | 12 | 19.9 s | 105.2 s | 157.0 s | 27.0 / 39.8 GB |

- **Grid against pages:** search 1.85x, build 3.3x (no global sort; 70% CPU instead
  of 57%), wall 1.9x, 12 GB less anonymous memory (the build's sort copies).
- **Shape:** R 4096 and m 10 are the best measured. Smaller blocks spread a window
  over more blocks; bigger ones scan more out-of-window peptides. Coarser cells scan
  longer runs; finer ones grow the directory (0.9 GB per dense slice at m 12, against
  0.22 GB at m 10).
- **Query order:** worth 5% with pages, nothing with the grid. The grid's lookups are
  independent loads (a cell start, then a short run), not dependent chains; why the
  order then doesn't matter was not measured further.

F9477 production, search phase:

| | pages | grid |
|---|---:|---:|
| run 10, warm page cache | 23.3 s | 20.4–21.3 s |
| run 11 | 125.0 s | 112.2 s |

Run 10 with pages measured 29–33 s while the page cache was cold (it streamed 9–15 GB
of spectra during the search); only the warm rerun is comparable. The monolithic
tryptic index is small (1,539 blocks) and run 11 spends most of its search in full
scoring, so the kernel's gain is diluted there. G23 tryptic search: 1.39 s → 0.91 s
monolithic, 3.6 s → 2.3 s sliced.

## Equal-count cells and 16-bit fragments (2026-10-16)

`c03cecb` left the dense class I slices bound by memory traffic at 16 threads: runs
of ~59 fragments (peaks fall where fragments are dense; equal-ppm cells held ~17 on
average), plus overscan, at 8 B per fragment. Two steps, each measured
(`plans/grid_compact_cells.md`):
- **a, equal-count cells** (cells from the fine-key histogram, kernel unchanged):
  scans fewer fragments where they are dense;
- **b, 16-bit residuals and block-local rows** (16-lane groups, the run rules above):
  4 B per fragment instead of 8.

Checks (`temp/fragment_grid/check_all.sh`, `cells/run.sh`, `compact/run.sh`,
necromerge2 root): unit tests for the cell table (monotone, covering, width cap,
`equal_count_cells` on a hand-made histogram), the residual round trip (every
generated (row, m/z bits) comes back from `fragments()`), and the kernel (dispatched,
scalar, AVX2) against the exact reference with windows inside one cell, across two and
many, and with edges exactly on fragment values, at R ∈ {1, 3, 8, 16, 64} ×
(f, T) ∈ {(7, 1), (10, 4), (14, 1), (14, 16), (16, 1000)}; `lanes_in`. Breaking the
middle-cell bounds, the position mask, the high-half mask bits, either scalar bound or
the AVX2 upper bound each fails them. Suites with the grid forced at five shapes
(including R 65,536). Data: small inputs, G23, F9477 runs 10/11 and every class I run
below, grid and pages, identical to `062a871`.

F9477 class I sliced, search seconds (and build), 16 threads unless noted:

| kernel | dense slices 1+2 | light slice 3 | search | build | wall | peak anon |
|---|---:|---:|---:|---:|---:|---:|
| `062a871` (grid) | 82.3 | 18.0 | 100.3 | 17.7 | 147.8 | 26.4 GB |
| `c03cecb` (branch-free, 8 lanes) | 79.3 | 16.4 | 95.7 | 19.1 | 146.3 | 26.4 GB |
| a, T 8 | 63.6 | 14.5 | 78.1 | 16.9 | 122.6 | 26.6 GB |
| a, T 16 | 66.1 | 16.2 | 82.3 | 20.2 | 131.8 | 26.4 GB |
| a, T 32 | 72.4 | 15.2 | 87.6 | 16.8 | 132.5 | 26.3 GB |
| **b, T 8** | **51.8** | 16.2 | **68.0** | 17.8 | **113.5** | **22.1 GB** |
| b, T 16 | 54.0 | 14.4 | 68.4 | 15.9 | 111.7 | 22.1 GB |
| b, f 16, T 4 | 54.5 | 14.2 | 68.7 | 16.8 | 112.5 | 22.6 GB |
| b, T 8, 8 threads | 62.9 | 20.2 | 83.1 | 19.6 | 136.5 | 22.2 GB |
| `c03cecb`, 8 threads | 81.5 | 24.5 | 106.0 | 16.9 | 158.0 | 26.1 GB |

- **Dense slices 1.59x faster than `062a871`**, search 1.47x, wall 1.30x; pages
  (`062a871`) searched these in 185.7 s.
- **Smaller cells help until the fine keys run out:** at T 8 the densest fine keys
  (f 14) already hold more than 8 fragments per block, so cells there are single fine
  keys (3,858 cells in slice 1, not 2 × 2,384). f 16 with T 4 doubles the directory
  (0.6–0.8 GB per dense slice) for no gain.
- **No longer purely bandwidth-bound:** 8 → 16 threads gives 1.21x on the dense slices
  again (`c03cecb`: 1.03x). The kernel now spends ~16 of its ~45 instructions per 16
  fragments on the position mask (`lanes_in`, branch-free `cmov`s); only a run's
  first and last group need it.
- **The light slice** (2–3M rows, ~480k queries) spends little time in the kernel;
  its 14–16 s varies with where the live memory budget puts the slice boundary.
- **Memory:** 4 B per fragment saves ~4 GB of the peak.

F9477 production search (monolithic tryptic): run 11 112.2 s (`062a871` grid) →
94.2 (`c03cecb`) → 88.1 (a, T 16) → 91.9 s (b, T 8); run 10 19.1 s (b) against
20.4 s (`062a871` grid). The tryptic index is small, so it gains least.

## Balance and class II (2026-10-16)

Each grid build also logs `fragment grid balance:` (`FragmentGrid::balance`): per tenth of
the blocks in mass order, the median fragments per block and the fragment-weighted cell
size (the size of the cell a fragment sits in, roughly what a window landing there
scans), plus the extremes.
- **Blocks:** fragments per block grow with mass, 1.2–1.4x within a class I slice; across
  class II's 16 slices from ~51k to ~150k (longer peptides), flat within each.
- **Cells:** the mean is ~8 fragments per block-cell (T 8), but the fragment-weighted size
  is 40–140, largest 1,834 (at R 1024, one cell held all 1,024 rows of its block), worse
  in light blocks. Likely clustered fragments: near-precursor ions (precursor minus one
  residue, ~20 values within a narrow-mass block) and short low-m/z ions with few
  distinct masses. Unconfirmed.

R re-swept with this kernel (class I wall): 4096 113.3 s, 2048 114.7 s, 1024 118.1 s.

**Class II** (9–25 aa, 620M rows, 16 slices): grid identical to pages (5,076,501 matched
fragments); wall 1,110 s → 380 s, peak anonymous 46.1 → 27.1 GB (the pages run also
re-read 47 GB from disk under memory pressure).

Search phases, class I, warm (in-source counters, temporary; 1,033 CPU-s): grid kernel
39%, spectrum read and processing 31% (65% in the light slice), full scoring 13%,
`add_matches` 8.5%, RT/IIM eviction and trim 6.8%, windows 2.7%, precursor binary search
1.8%. A cheaper run-edge mask (each run's first and last mask computed once, chosen by
`cmov`: ~38 instead of ~45 instructions per 16 fragments) changed class I search by 2% and
class II not at all, so it was not kept: the kernel is no longer bound by instructions.

## Predicted RT as a third axis (2026-10-16)

With `--predicted-rt` (and `rt_tol_sec`) and the grid layout, each block's fragments are
also split by their peptide's predicted RT:

```
block (rows in mass order) -> RT slab -> m/z cell -> entries
directory per block: slabs * cells + 1 starts; slab s, cell c at entry s * cells + c
```

- **Slabs:** `RtSlabbing { predicted, slabs }` (`fragment_grid.rs`) gives one predicted RT
  per row, the same array `Scorer::predicted_rt` reads. Slab bounds are equal-count
  quantiles of a sample of the finite predictions (`rt_slab_bounds`); a row with
  predicted RT `rt` is in slab `rt_bounds.partition_point(|&b| b <= rt)`. Rows without a
  prediction, or with a non-finite one, go to one extra slab that every spectrum reads.
- **How many:** `RtSlabbing::derived_slabs`: the predictions' span over four widths of a
  spectrum's compatible window at mid-span (`Runner::rt_slabbing`), 1 to 64;
  `SAGE_GRID_RT_SLABS` overrides it. Cells then target T fragments per (block, slab, cell)
  (`cell_fragments * blocks * slabs` per cell), so the directory stays about as large.
- **Where it is built:** sliced search, per slice from its rows' predictions
  (`build_from_peptides_with_rt`); monolithic search, by rebuilding the grid once the
  predictions are resolved (`Runner::run`; the database is built before predictions are
  read).
- **Query:** `Scorer::compatible_predicted_rt` (the range eviction keeps, now shared with
  the kernel) → `FragmentGrid::slab_visits`: the RT slabs overlapping it and the
  no-prediction slab. Rows in skipped slabs have predicted RTs outside the range, so
  sage's eviction would have removed them; rows of the visited slabs that fall outside are
  counted and then evicted as before.
- **Exact without changing sage:** a window whose every match is evicted leaves sage with
  an all-empty preliminary array trimmed to ≤ 50 slots, while pruning leaves it empty and
  untrimmed. Pooled trims see the same first k slots either way and empty slots never
  enter `bounded_min_heapify`'s heap, so the pooled results, and the output, are the same.
- **Tests:** `rt_slabs_keep_every_match_sage_keeps_and_add_none`: random rows with
  predictions, none, NaN and infinity; RT windows inside, around, outside, inverted, NaN,
  and on every slab bound; R ∈ {1, 7, 64} × slabs ∈ {1, 2, 5, 13}: the kernel's matches
  contain every match sage's RT check keeps and nothing the unpruned search lacks; scalar
  equals dispatched. Dropping the last visited RT slab, skipping the no-prediction slab or
  a strict bound in the build each fails it.

RT slabs on F9477 run 11 (monolithic tryptic, `--predicted-rt`, R 4096 unless noted),
search seconds, every run identical to pages:

| slabs (RT + no prediction) | cells | search |
|---|---:|---:|
| no RT axis | 10,146 | 75.9 |
| 1 + 1 | 10,146 | 78.2 |
| **4 + 1** | 3,513 | **69.9** |
| 17 + 1 (span over one window) | 1,058 | 71.5 |
| 32 + 1 | 690 | 76.0 |
| 17 + 1, R 256 | 588 | 70.5 (directory 1.04 GB) |

The kernel is about a fifth of this search, so the RT axis gains ~8%. More slabs leave
fewer fragments per (block, slab) and so fewer, wider cells; past ~4 slabs the m/z
waste outweighs the rows pruned, hence `derived_slabs` divides the span by four windows
(4 slabs here). One dataset only.

