# Fragment grid: precursor blocks × fragment cells (2026-10-16)

Necromerge2 plan: `plans/fragment_grid.md`. A second layout of the fragment index
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
cell_start[block][c] -> offset (relative to the block) into the (mz, peptide) columns
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
- **Storage:** the same `fragment_mzs`/`fragment_peptide_indices` columns, in grid
  order, `min_value` empty, plus `IndexedDatabase::grid` (block starts and per-block
  `cells + 1` cell starts, `u32` relative to the block).

## Build

`Parameters::build_from_peptides` with `fragment_layout = Grid` calls `build_grid`.
The fragments come from `Parameters::indexed_fragments`, shared with the page layout,
so their m/z bits are identical. `FragmentGrid::build` never holds more than the
columns and the directory:
1. **Count:** per block, in parallel, generate its fragments into thread-local
   scratch; keep the count and the lowest/highest cell key.
2. **Offsets:** block starts by prefix sum; the global key range gives `first_cell`
   and the cell count.
3. **Place:** per block, in parallel, into its own disjoint range of the columns:
   generate again, count per cell, prefix sum (the block's cell starts), place in
   generation order (a stable counting sort).

No global sort. The page layout instead sorts every fragment by m/z, then sorts each
page by peptide, and peaks at about 16 B per fragment for those copies.

## Search

`IndexedQuery::page_search_batch` hands the grid the same inputs the pages get: the
exact candidate rows `[exact_pre_idx_lo, pre_idx_hi)` and the fragment windows.
`FragmentGrid::search`:
- per window, once: its `(lo, hi)` and the cells holding it;
- per block overlapping the candidate rows (outer loop), per window (inner loop):
  scan the run `cell_start[b][first] .. cell_start[b][last + 1]` with the AVX2/scalar
  `lo <= mz <= hi` comparison (`scan_masses`); in a block not wholly inside the
  candidate rows, also check the peptide id.

**Exact by construction.** The cells partition each block and the key is monotone, so
every fragment in a window lies in the scanned cells, and the comparison and the id
check select exactly the (window, fragment) pairs the pages select. Only the order
they are reported in differs. `matched_peaks_with_isotope` only counts per candidate,
so the preliminary array, RT/IIM eviction and `trim_hits` see the same values.
`page_search` (one window; only an integration test uses it) goes through the same
kernel.

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
