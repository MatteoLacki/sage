# Fragment grid: precursor blocks × fragment cells (2026-10-16)

Necromerge2 plans: `plans/fragment_grid.md`, `plans/grid_match_simd.md` (the
branch-free match kernel, below). A second layout of the fragment index
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
                 fragment cells c = (mz.to_bits() >> (23 - m)) - first_cell
                 c0     c1     c2    ...
block 0  rows 0..R    [run]  [run]  [run]  ...
block 1  rows R..2R   [run]  [run]  [run]  ...
cell_start[block][c] -> fragment offset (relative to the block) into its 8-fragment groups
```

- **Blocks:** R consecutive peptide rows. Rows are mass-sorted, so a block is a
  precursor-mass band; equal row counts keep blocks equally full where peptide
  density varies.
- **Cells:** for non-negative floats the bit pattern is monotone in the value, and its
  top bits are the exponent and leading mantissa bits. Keeping m mantissa bits splits
  each octave into 2^m equal cells, about constant ppm wide (m = 10: 488–977 ppm).
  One shift, no logarithm. `first_cell` is the lightest fragment's key; window keys
  are clamped to the grid, and a window whose high end is not positive holds nothing.
- **Inside a cell:** generation order, i.e. peptide order. Nothing relies on it.
- **Storage:** `FragmentGrid` owns the fragments as `FragmentGroup`s, 64 B each (one
  cache line, `align(64)`): 8 m/z (`f32`), then their 8 peptide rows (`u32`). Each
  block starts on a group and is padded to whole groups with `FragmentGroup::PAD`
  (NaN m/z, which no window matches). Also block group starts and per-block
  `cells + 1` cell starts (`u32` fragment offsets relative to the block).
  `IndexedDatabase::grid` holds it; `fragment_mzs`/`fragment_peptide_indices` and
  `min_value` stay empty. `size()` and `fragments()` read the grid, pads excluded
  (`database_size` and `dump_fragment_index` read the page columns, so they see
  nothing for a grid).

## Build

`Parameters::build_from_peptides` with `fragment_layout = Grid` calls `build_grid`.
The fragments come from `Parameters::indexed_fragments`, shared with the page layout,
so their m/z bits are identical. `FragmentGrid::build` never holds more than the
columns and the directory:
1. **Count:** per block, in parallel, generate its fragments into thread-local
   scratch; keep the count and the lowest/highest cell key.
2. **Offsets:** block starts by prefix sum; the global key range gives `first_cell`
   and the cell count.
3. **Place:** the groups are allocated filled with pads (a parallel collect); per
   block, in parallel, into its own disjoint groups: generate again, count per cell,
   prefix sum (the block's cell starts), place in generation order (a stable counting
   sort) into lane `offset % 8` of group `offset / 8`.

No global sort. The page layout instead sorts every fragment by m/z, then sorts each
page by peptide, and peaks at about 16 B per fragment for those copies.

## Search

`matched_peaks_with_isotope` asks `IndexedQuery::matched_candidates` for the matched
candidates' offsets from `candidate_rows().start` (the exact rows
`[exact_pre_idx_lo, pre_idx_hi)`), one entry per (window, fragment) match, and counts
them into `preliminary` (`InitialHits::add_matches`). The grid answers with
`FragmentGrid::matched_candidates`; pages run their own `page_search_batch` with a
closure that pushes offsets. Per call:
- per window, once: its `(lo, hi)` and the cells holding it;
- per block overlapping the candidate rows (outer loop), per window (inner loop): the
  run of groups covering `cell_start[b][first] .. cell_start[b][last + 1]`, whole
  groups (below).

**Exact by construction.** The cells partition each block and the key is monotone, so
every fragment in a window lies in the scanned cells. A run's first and last groups
can also hold fragments of neighbouring cells (keys outside the window's, so m/z
outside it) and pads (NaN), never another block's. The m/z comparison and the
candidate check therefore select exactly the (window, fragment) pairs the pages
select; only their order differs, and counting doesn't depend on it.

`FragmentGrid::search` (the per-fragment callback, reached through
`page_search_batch` and `page_search`) visits each run's fragments one at a time,
exactly; the tests use it as the reference.

## Branch-free match handling (2026-10-16)

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
| `SAGE_GRID_ROWS` | R | 4096 |
| `SAGE_GRID_CELL_BITS` | m, 0–23 | 10 |
| `SAGE_QUERY_ORDER` | `mass`, `file` (sliced search) | `mass` |

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
