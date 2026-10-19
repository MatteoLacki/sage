//! Precursor-major fragment index: blocks of consecutive peptide rows, each split into
//! fragment m/z cells (necromerge2 `plans/fragment_grid.md`, `docs/ai/fragment_grid.md`).

use crate::database::{PeptideIx, Theoretical};
use crate::mass::Tolerance;
use crate::predicted::PredictedValues;
use rayon::prelude::*;
use std::cell::RefCell;
use std::ops::Range;
use std::sync::atomic::{AtomicBool, Ordering as AtomicOrdering};

/// How `Parameters::build_from_peptides` lays out the fragment index.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum FragmentLayout {
    /// Fragment-mass pages, each sorted by peptide (sage's own).
    #[default]
    Pages,
    Grid(GridShape),
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GridShape {
    pub rows_per_block: usize,
    /// Mantissa bits kept in a fine key: 2^fine_bits fine keys per octave of m/z. Cells
    /// are runs of consecutive fine keys.
    pub fine_bits: u32,
    /// Fragments a cell holds per block, on average over the grid.
    pub cell_fragments: usize,
}

impl Default for GridShape {
    fn default() -> Self {
        GridShape {
            rows_per_block: 4096,
            fine_bits: 14,
            cell_fragments: 8,
        }
    }
}

/// Coarser fine keys would span more than 2^16 bit patterns each, the most one cell may.
const MIN_FINE_BITS: u32 = 7;
/// Bounds the census histogram (2^fine_bits bins per octave of `FINE_KEY_MZ`).
const MAX_FINE_BITS: u32 = 16;
/// The m/z range a fragment must lie in to be placed (the census histogram's range).
const FINE_KEY_MZ: Range<f32> = 16.0..32768.0;

impl FragmentLayout {
    /// `SAGE_FRAGMENT_LAYOUT=pages|grid`, with `SAGE_GRID_ROWS`, `SAGE_GRID_FINE_BITS`
    /// and `SAGE_GRID_CELL_FRAGMENTS` shaping the grid.
    pub fn from_environment() -> Self {
        match std::env::var("SAGE_FRAGMENT_LAYOUT").as_deref() {
            Err(_) | Ok("pages") => FragmentLayout::Pages,
            Ok("grid") => {
                assert!(
                    std::env::var("SAGE_GRID_CELL_BITS").is_err(),
                    "SAGE_GRID_CELL_BITS (equal-ppm cells) is replaced by SAGE_GRID_FINE_BITS \
                     and SAGE_GRID_CELL_FRAGMENTS (equal-count cells)"
                );
                let default = GridShape::default();
                let shape = GridShape {
                    rows_per_block: environment_number("SAGE_GRID_ROWS")
                        .unwrap_or(default.rows_per_block),
                    fine_bits: environment_number("SAGE_GRID_FINE_BITS")
                        .map(|bits| bits as u32)
                        .unwrap_or(default.fine_bits),
                    cell_fragments: environment_number("SAGE_GRID_CELL_FRAGMENTS")
                        .unwrap_or(default.cell_fragments),
                };
                assert!(
                    shape.rows_per_block >= 1,
                    "SAGE_GRID_ROWS must be at least 1"
                );
                assert!(
                    (MIN_FINE_BITS..=MAX_FINE_BITS).contains(&shape.fine_bits),
                    "SAGE_GRID_FINE_BITS must be {MIN_FINE_BITS} to {MAX_FINE_BITS}"
                );
                assert!(
                    shape.cell_fragments >= 1,
                    "SAGE_GRID_CELL_FRAGMENTS must be at least 1"
                );
                FragmentLayout::Grid(shape)
            }
            Ok(other) => panic!("SAGE_FRAGMENT_LAYOUT must be pages or grid, not {other}"),
        }
    }
}

/// Predicted retention times by peptide row (per-run calibrated, in the units of the
/// spectra's retention times), and how many RT slabs to split each block into.
#[derive(Clone, Copy, Debug)]
pub struct RtSlabbing<'a> {
    pub predicted: PredictedValues<'a>,
    pub slabs: usize,
}

impl RtSlabbing<'_> {
    /// `SAGE_GRID_RT_SLABS`, else the span of `predicted` over four times `window`, the
    /// width of a spectrum's compatible window; 1 to 64. More slabs prune more rows but
    /// leave fewer fragments per (block, slab) and so wider cells; on F9477 run 11 the
    /// span over one window (17 slabs) was slower than over four (4 slabs).
    pub fn derived_slabs(predicted: PredictedValues, window: f32) -> usize {
        if let Some(slabs) = environment_number("SAGE_GRID_RT_SLABS") {
            assert!(slabs >= 1, "SAGE_GRID_RT_SLABS must be at least 1");
            return slabs;
        }
        let finite = (0..predicted.len())
            .map(|row| predicted.get(row))
            .filter(|rt| rt.is_finite());
        let (low, high) = finite.fold((f32::INFINITY, f32::NEG_INFINITY), |(low, high), rt| {
            (low.min(rt), high.max(rt))
        });
        match high > low && window > 0.0 {
            true => ((high - low) / (4.0 * window)).floor().clamp(1.0, 64.0) as usize,
            false => 1,
        }
    }
}

fn environment_number(name: &str) -> Option<usize> {
    let value = std::env::var(name).ok()?;
    Some(
        value
            .parse()
            .unwrap_or_else(|_| panic!("{name} must be a whole number, not {value}")),
    )
}

const LANES: usize = 16;
/// Lanes per half group: one `LEFT_PACK` entry each.
const HALF: usize = 8;

/// Sixteen fragments in one cache line: each one's m/z as its offset from its cell's
/// first bit pattern (`m/z bits = cell_first_bits + residual`), then its peptide row
/// within the block. Lanes past a block's last fragment are never read.
#[repr(C, align(64))]
#[derive(Clone, Copy, Debug)]
pub struct FragmentGroup {
    residuals: [u16; LANES],
    rows: [u16; LANES],
}

impl FragmentGroup {
    const EMPTY: FragmentGroup = FragmentGroup {
        residuals: [u16::MAX; LANES],
        rows: [u16::MAX; LANES],
    };
}

/// For each 8-bit lane mask, the set lanes in order: the `permutevar8x32` indices
/// that move a half group's matching lanes to the front.
static LEFT_PACK: [[u32; HALF]; 256] = left_pack_table();

const fn left_pack_table() -> [[u32; HALF]; 256] {
    let mut table = [[0u32; HALF]; 256];
    let mut mask = 0;
    while mask < 256 {
        let mut lane = 0;
        let mut packed = 0;
        while lane < HALF {
            if mask & (1 << lane) != 0 {
                table[mask][packed] = lane as u32;
                packed += 1;
            }
            lane += 1;
        }
        mask += 1;
    }
    table
}

/// The fragments, in groups of 16, and where each (block, cell) run of them starts.
#[derive(Debug, Default)]
pub struct FragmentGrid {
    rows: usize,
    rows_per_block: usize,
    cell_fragments: usize,
    /// `mz.to_bits() >> shift` is a fragment's fine key.
    shift: u32,
    /// The lightest fragment's fine key: `cell_of_fine[0]`.
    first_fine: u32,
    /// The cell of every fine key from `first_fine` to the heaviest fragment's.
    cell_of_fine: Vec<u16>,
    /// Each cell's lowest fine key.
    cell_first_fine: Vec<u32>,
    cells: usize,
    /// Slabs per block: `rt_bounds.len() + 1` RT slabs and one for rows without a usable
    /// prediction, or a single slab holding everything when there is no RT axis.
    slabs: usize,
    /// Upper bounds of the RT slabs but the last, ascending: a row with predicted RT `rt`
    /// is in slab `rt_bounds.partition_point(|&bound| bound <= rt)`.
    rt_bounds: Vec<f32>,
    rt_axis: bool,
    /// Fragments, pads excluded.
    fragments: usize,
    /// Each block's first group, plus the total.
    block_groups: Vec<usize>,
    /// Per block, `slabs * cells + 1` fragment offsets relative to the block's first
    /// group: slab `s`, cell `c` starts at entry `s * cells + c`.
    cell_starts: Vec<u32>,
    groups: Vec<FragmentGroup>,
}

/// One searched window: its fragment bounds, the cells holding them, and the bounds as
/// residuals in its first and last cell.
#[derive(Clone, Copy, Debug)]
struct CellWindow {
    lo: f32,
    hi: f32,
    first: usize,
    last: usize,
    residual_lo: u16,
    residual_hi: u16,
}

impl CellWindow {
    /// The window's fragments in a block as three runs (some empty), each with the
    /// residual bounds its fragments must meet: the first cell (`>= residual_lo`), the
    /// cells between (inside the window by construction), the last cell
    /// (`<= residual_hi`); one run with both bounds when the cells coincide.
    fn runs(&self, starts: &[u32]) -> [(Range<usize>, u16, u16); 3] {
        let start = |cell: usize| starts[cell] as usize;
        if self.first == self.last {
            let run = start(self.first)..start(self.first + 1);
            return [
                (run, self.residual_lo, self.residual_hi),
                (0..0, 0, 0),
                (0..0, 0, 0),
            ];
        }
        [
            (
                start(self.first)..start(self.first + 1),
                self.residual_lo,
                u16::MAX,
            ),
            (start(self.first + 1)..start(self.last), 0, u16::MAX),
            (start(self.last)..start(self.last + 1), 0, self.residual_hi),
        ]
    }
}

/// Per-thread scratch for building and searching.
#[derive(Default)]
struct GridScratch {
    fragments: Vec<Theoretical>,
    cursors: Vec<u32>,
    /// The RT slab of each row of the block being built.
    row_slabs: Vec<u16>,
    windows: Vec<CellWindow>,
}

thread_local! {
    static GRID_SCRATCH: RefCell<GridScratch> = RefCell::new(GridScratch::default());
}

/// Blocks one sampled block stands for when a grid's cells are chosen.
const CENSUS_EVERY: usize = 4;

/// `SAGE_GRID_CENSUS_EVERY` (1: count every block), else [`CENSUS_EVERY`].
pub fn census_every() -> usize {
    environment_number("SAGE_GRID_CENSUS_EVERY")
        .unwrap_or(CENSUS_EVERY)
        .max(1)
}

impl FragmentGrid {
    /// Lay out the fragments of peptide rows `0..rows`, generated by `fragments_of(row,
    /// out)` in a fixed order, as a grid. Each block is generated twice (count, then
    /// place), so nothing beyond the grid is held.
    pub fn build(
        rows: usize,
        shape: GridShape,
        rt: Option<RtSlabbing>,
        fragments_of: impl Fn(usize, &mut Vec<Theoretical>) + Sync,
    ) -> FragmentGrid {
        Self::build_sampled(rows, shape, rt, &fragments_of, None::<&fn(usize) -> usize>, 1)
    }

    /// `build`, with the cells chosen from the fragments of every `every`-th block only
    /// (and the last), and block sizes from `fragment_count(row)`: a block is generated
    /// once more only when sampled. Any increasing cell table finds the same matches,
    /// so this changes how evenly cells fill, never the result. A fragment outside the
    /// sampled range (widened by a margin) falls back to counting every block.
    pub fn build_sampled<F, C>(
        rows: usize,
        shape: GridShape,
        rt: Option<RtSlabbing>,
        fragments_of: &F,
        fragment_count: Option<&C>,
        every: usize,
    ) -> FragmentGrid
    where
        F: Fn(usize, &mut Vec<Theoretical>) + Sync,
        C: Fn(usize) -> usize + Sync,
    {
        let every = match fragment_count {
            Some(_) => every.max(1),
            None => 1,
        };
        if let Some(rt) = &rt {
            assert_eq!(rt.predicted.len(), rows, "one predicted RT per row");
        }
        let rt_bounds = rt.as_ref().map(rt_slab_bounds).unwrap_or_default();
        let rt_slabs = rt.map_or(0, |_| rt_bounds.len() + 1);
        let slabs = rt_slabs + 1;
        let slab_of_row = |row: usize| -> usize {
            let Some(rt) = rt else {
                return 0;
            };
            let predicted = rt.predicted.get(row);
            match predicted.is_finite() {
                true => rt_bounds.partition_point(|&bound| bound <= predicted),
                false => rt_slabs,
            }
        };
        let slabs = if rt.is_some() { slabs } else { 1 };
        assert!(rows <= u32::MAX as usize, "peptide rows must fit u32");
        assert!(
            shape.rows_per_block <= 1 << 16,
            "rows within a block must fit u16"
        );
        let shift = 23 - shape.fine_bits;
        let histogram_base = FINE_KEY_MZ.start.to_bits() >> shift;
        let histogram_bins = ((FINE_KEY_MZ.end.to_bits() >> shift) - histogram_base) as usize;
        let blocks = (rows + shape.rows_per_block - 1) / shape.rows_per_block;
        let block_rows = |block: usize| {
            block * shape.rows_per_block..((block + 1) * shape.rows_per_block).min(rows)
        };
        let generate = |block: usize, out: &mut Vec<Theoretical>| {
            out.clear();
            for row in block_rows(block) {
                fragments_of(row, out);
            }
        };

        let splits = 4 * rayon::current_num_threads();
        let sampled = |block: usize| block % every == 0 || block + 1 == blocks;
        let (mut census, histogram) = (0..blocks)
            .into_par_iter()
            .with_min_len((blocks / splits).max(1))
            .filter(|&block| sampled(block))
            .fold(
                || (Vec::new(), vec![0u32; histogram_bins]),
                |(mut census, mut histogram), block| {
                    GRID_SCRATCH.with(|scratch| {
                        let fragments = &mut scratch.borrow_mut().fragments;
                        generate(block, fragments);
                        let mut lowest = u32::MAX;
                        let mut highest = 0;
                        for fragment in fragments.iter() {
                            assert!(
                                FINE_KEY_MZ.contains(&fragment.fragment_mz),
                                "fragment m/z {} is outside the grid's {FINE_KEY_MZ:?}",
                                fragment.fragment_mz
                            );
                            let key = fragment.fragment_mz.to_bits() >> shift;
                            lowest = lowest.min(key);
                            highest = highest.max(key);
                            histogram[(key - histogram_base) as usize] += 1;
                        }
                        census.push((block, fragments.len(), lowest, highest));
                    });
                    (census, histogram)
                },
            )
            .reduce(
                || (Vec::new(), vec![0u32; histogram_bins]),
                |(mut census, mut histogram), (other_census, other_histogram)| {
                    census.extend(other_census);
                    for (count, other) in histogram.iter_mut().zip(other_histogram) {
                        *count += other;
                    }
                    (census, histogram)
                },
            );
        census.sort_unstable_by_key(|&(block, ..)| block);
        let block_counts: Vec<usize> = match fragment_count {
            Some(fragment_count) if every > 1 => (0..blocks)
                .into_par_iter()
                .map(|block| block_rows(block).map(fragment_count).sum())
                .collect(),
            _ => census.iter().map(|&(_, count, ..)| count).collect(),
        };

        let lowest = census
            .iter()
            .map(|&(_, _, lowest, _)| lowest)
            .min()
            .unwrap_or(0);
        let highest = census
            .iter()
            .map(|&(_, _, _, highest)| highest)
            .max()
            .unwrap_or(0);
        let (first_fine, last_fine) = match every {
            1 => (lowest, highest),
            _ => {
                let margin = 1u32 << (shape.fine_bits - 5);
                let histogram_last = histogram_base + histogram_bins as u32 - 1;
                (
                    lowest.saturating_sub(margin).max(histogram_base),
                    (highest + margin).min(histogram_last),
                )
            }
        };
        let sampled_fragments: u64 = census.iter().map(|&(_, count, ..)| count as u64).sum();
        let all_fragments: u64 = block_counts.iter().map(|&count| count as u64).sum();
        let target = (shape.cell_fragments * blocks * rt_slabs.max(1)) as u64;
        let (cell_of_fine, cell_first_fine) = match sampled_fragments > 0 {
            true => equal_count_cells(
                &histogram[(first_fine - histogram_base) as usize
                    ..=(last_fine - histogram_base) as usize],
                first_fine,
                (target * sampled_fragments / all_fragments).max(1),
                1 << (shape.fine_bits - MIN_FINE_BITS),
            ),
            false => (Vec::new(), Vec::new()),
        };
        let cells = cell_first_fine.len();
        let mut block_groups = Vec::with_capacity(blocks + 1);
        block_groups.push(0);
        for &count in &block_counts {
            assert!(
                count <= u32::MAX as usize,
                "a grid block holds too many fragments"
            );
            block_groups.push(block_groups.last().unwrap() + (count + LANES - 1) / LANES);
        }
        let fragments = all_fragments as usize;

        let mut groups: Vec<FragmentGroup> = (0..*block_groups.last().unwrap())
            .into_par_iter()
            .map(|_| FragmentGroup::EMPTY)
            .collect();
        let directory = slabs * cells + 1;
        let mut cell_starts = vec![0u32; blocks * directory];
        let out_of_range = AtomicBool::new(false);
        split_by_starts(&mut groups, &block_groups)
            .into_par_iter()
            .zip(cell_starts.par_chunks_mut(directory))
            .zip(block_counts.par_iter())
            .enumerate()
            .for_each(|(block, ((block_groups, starts), &count))| {
                GRID_SCRATCH.with(|scratch| {
                    let GridScratch {
                        fragments,
                        cursors,
                        row_slabs,
                        ..
                    } = &mut *scratch.borrow_mut();
                    generate(block, fragments);
                    assert_eq!(fragments.len(), count, "block {block} changed size");
                    let outside_the_cells = |fragment: &Theoretical| {
                        let key = fragment.fragment_mz.to_bits() >> shift;
                        !FINE_KEY_MZ.contains(&fragment.fragment_mz)
                            || key < first_fine
                            || key > last_fine
                    };
                    if every > 1 && fragments.iter().any(outside_the_cells) {
                        out_of_range.store(true, AtomicOrdering::Relaxed);
                        return;
                    }
                    let first_row = block * shape.rows_per_block;
                    row_slabs.clear();
                    row_slabs.extend(block_rows(block).map(|row| slab_of_row(row) as u16));
                    let cell_of =
                        |bits: u32| cell_of_fine[((bits >> shift) - first_fine) as usize] as usize;
                    let key_of = |fragment: &Theoretical| {
                        let slab =
                            row_slabs[fragment.peptide_index.0 as usize - first_row] as usize;
                        slab * cells + cell_of(fragment.fragment_mz.to_bits())
                    };
                    for fragment in fragments.iter() {
                        starts[key_of(fragment) + 1] += 1;
                    }
                    for cell in 1..starts.len() {
                        starts[cell] += starts[cell - 1];
                    }
                    cursors.clear();
                    cursors.extend_from_slice(starts);
                    for fragment in fragments.iter() {
                        let bits = fragment.fragment_mz.to_bits();
                        let cell = cell_of(bits);
                        let cursor = &mut cursors[key_of(fragment)];
                        let group = &mut block_groups[*cursor as usize / LANES];
                        let lane = *cursor as usize % LANES;
                        let residual = bits - (cell_first_fine[cell] << shift);
                        debug_assert!(residual <= u16::MAX as u32, "cell wider than 2^16");
                        group.residuals[lane] = residual as u16;
                        group.rows[lane] = (fragment.peptide_index.0 as usize - first_row) as u16;
                        *cursor += 1;
                    }
                });
            });
        if out_of_range.load(AtomicOrdering::Relaxed) {
            log::info!(
                "fragment grid: a fragment lies outside the sampled census's range, counting every block"
            );
            return Self::build_sampled(rows, shape, rt, fragments_of, fragment_count, 1);
        }

        FragmentGrid {
            rows,
            rows_per_block: shape.rows_per_block,
            cell_fragments: shape.cell_fragments,
            shift,
            first_fine,
            cell_of_fine,
            cell_first_fine,
            cells,
            slabs,
            rt_bounds,
            rt_axis: rt.is_some(),
            fragments,
            block_groups,
            cell_starts,
            groups,
        }
    }

    pub fn describe(&self) -> String {
        format!(
            "{} blocks of {} rows x {} RT slabs x {} cells (~{} fragments each, {}-bit fine keys), {} fragments in {} groups, directory {:.1} MB",
            self.block_groups.len() - 1,
            self.rows_per_block,
            self.slabs,
            self.cells,
            self.cell_fragments,
            23 - self.shift,
            self.fragments,
            self.groups.len(),
            (self.cell_starts.len() * std::mem::size_of::<u32>()) as f64 / 1e6
        )
    }

    /// How evenly the fragments spread, in mass order: per tenth of the blocks, the median
    /// fragments per block and the fragment-weighted mean cell size (the size of the cell
    /// a fragment sits in: what a window landing on it scans), plus the overall extremes.
    pub fn balance(&self) -> String {
        let blocks = self.block_groups.len().saturating_sub(1);
        if blocks == 0 || self.cells == 0 {
            return "empty".into();
        }
        let directory = self.slabs * self.cells + 1;
        let block_fragments = |block: usize| {
            let starts = &self.cell_starts[block * directory..][..directory];
            starts[directory - 1] as usize
        };
        let tenths = blocks.min(10);
        let mut per_block = Vec::with_capacity(tenths);
        let mut weighted_cell = Vec::with_capacity(tenths);
        for tenth in 0..tenths {
            let range = tenth * blocks / tenths..(tenth + 1) * blocks / tenths;
            let mut counts: Vec<usize> = range.clone().map(block_fragments).collect();
            counts.sort_unstable();
            per_block.push(counts[counts.len() / 2]);
            let (mut squares, mut total) = (0u64, 0u64);
            for block in range {
                let starts = &self.cell_starts[block * directory..][..directory];
                for pair in starts.windows(2) {
                    let size = u64::from(pair[1] - pair[0]);
                    squares += size * size;
                    total += size;
                }
            }
            weighted_cell.push(squares as f64 / total.max(1) as f64);
        }
        let all: Vec<usize> = (0..blocks).map(block_fragments).collect();
        let largest_cell = (0..blocks)
            .flat_map(|block| {
                let starts = &self.cell_starts[block * directory..][..directory];
                starts.windows(2).map(|pair| pair[1] - pair[0]).max()
            })
            .max()
            .unwrap_or(0);
        format!(
            "fragments per block {}..{} (median by tenth of blocks {:?}); fragment-weighted cell size by tenth {:?}, largest cell {}",
            all.iter().min().unwrap(),
            all.iter().max().unwrap(),
            per_block,
            weighted_cell
                .iter()
                .map(|size| (size * 10.0).round() / 10.0)
                .collect::<Vec<_>>(),
            largest_cell
        )
    }

    /// Fragments, pads excluded.
    pub fn len(&self) -> usize {
        self.fragments
    }

    pub fn is_empty(&self) -> bool {
        self.fragments == 0
    }

    /// Every fragment, in grid order.
    pub fn fragments(&self) -> impl Iterator<Item = Theoretical> + '_ {
        (0..self.block_groups.len().saturating_sub(1)).flat_map(move |block| {
            let (groups, directory) = self.block(block);
            (0..self.slabs).flat_map(move |slab| {
                let starts = self.slab_starts(directory, slab);
                (0..self.cells).flat_map(move |cell| {
                    (starts[cell] as usize..starts[cell + 1] as usize)
                        .map(move |fragment| self.fragment(block, groups, cell, fragment))
                })
            })
        })
    }

    /// The fragment at `fragment` of block `block` (whose groups are `groups`), in cell
    /// `cell`: its exact m/z and table row.
    fn fragment(
        &self,
        block: usize,
        groups: &[FragmentGroup],
        cell: usize,
        fragment: usize,
    ) -> Theoretical {
        let group = &groups[fragment / LANES];
        let lane = fragment % LANES;
        Theoretical {
            peptide_index: PeptideIx(
                (block * self.rows_per_block + group.rows[lane] as usize) as u32,
            ),
            fragment_mz: f32::from_bits(self.cell_first_bits(cell) + group.residuals[lane] as u32),
        }
    }

    fn cell_first_bits(&self, cell: usize) -> u32 {
        self.cell_first_fine[cell] << self.shift
    }

    /// Block `block`'s candidate rows, local to the block, and what turns a local row
    /// into an offset from `candidates.start` (added with wrapping).
    fn local_candidates(&self, block: usize, candidates: &Range<usize>) -> (Range<u32>, u32) {
        let first_row = block * self.rows_per_block;
        let end_row = (first_row + self.rows_per_block).min(self.rows);
        let local = (candidates.start.max(first_row) - first_row) as u32
            ..(candidates.end.min(end_row) - first_row) as u32;
        (
            local,
            (first_row as u32).wrapping_sub(candidates.start as u32),
        )
    }

    fn blocks_holding(&self, candidates: &Range<usize>) -> std::ops::RangeInclusive<usize> {
        candidates.start / self.rows_per_block..=(candidates.end - 1) / self.rows_per_block
    }

    /// Block `block`'s groups and its whole directory (`slabs * cells + 1` starts).
    fn block(&self, block: usize) -> (&[FragmentGroup], &[u32]) {
        let groups = &self.groups[self.block_groups[block]..self.block_groups[block + 1]];
        let directory = self.slabs * self.cells + 1;
        let starts = &self.cell_starts[block * directory..][..directory];
        (groups, starts)
    }

    /// The `cells + 1` cell starts of slab `slab` within a block's directory.
    fn slab_starts<'a>(&self, directory: &'a [u32], slab: usize) -> &'a [u32] {
        &directory[slab * self.cells..][..self.cells + 1]
    }

    /// The slabs a spectrum whose compatible predicted RTs are `rt` (`None`: unknown)
    /// must read, as two ranges: the RT slabs overlapping the range, and the slab of rows
    /// without a usable prediction (sage keeps those candidates). Every slab when there
    /// is no RT axis or no range.
    fn slab_visits(&self, rt: Option<(f32, f32)>) -> [Range<usize>; 2] {
        let (Some((low, high)), true) = (rt, self.rt_axis) else {
            return [0..self.slabs, 0..0];
        };
        let no_prediction = self.slabs - 1;
        let some_rt_is_compatible = low <= high;
        if !some_rt_is_compatible {
            return [0..0, no_prediction..self.slabs];
        }
        let first = self.rt_bounds.partition_point(|&bound| bound <= low);
        let last = self.rt_bounds.partition_point(|&bound| bound <= high);
        [first..last + 1, no_prediction..self.slabs]
    }

    /// Each window's fragment bounds and cells, skipping windows outside the grid. A
    /// window's low end lies in its first cell or below it (then every residual
    /// qualifies), its high end in its last cell or above it; residuals of a cell are
    /// the bit patterns' order, so the bounds are exact.
    fn cell_windows(&self, windows: &[(f32, Tolerance)], out: &mut Vec<CellWindow>) {
        out.clear();
        for &(mass, tolerance) in windows {
            let (lo, hi) = tolerance.theoretical_range(mass);
            let Some((first, last)) = self.cells_between(lo, hi) else {
                continue;
            };
            let lo_bits = if lo > 0.0 { lo.to_bits() } else { 0 };
            let residual_lo = lo_bits.saturating_sub(self.cell_first_bits(first));
            let residual_hi = hi.to_bits().saturating_sub(self.cell_first_bits(last));
            debug_assert!(residual_lo <= u16::MAX as u32);
            out.push(CellWindow {
                lo,
                hi,
                first,
                last,
                residual_lo: residual_lo as u16,
                residual_hi: residual_hi.min(u16::MAX as u32) as u16,
            });
        }
    }

    /// The cells holding every fragment with `lo <= m/z <= hi`; `None` when none can.
    fn cells_between(&self, lo: f32, hi: f32) -> Option<(usize, usize)> {
        let window_can_hold_a_fragment = hi >= lo && hi > 0.0;
        if self.cells == 0 || !window_can_hold_a_fragment {
            return None;
        }
        let last_fine = self.first_fine + (self.cell_of_fine.len() - 1) as u32;
        let lo_key = match lo > 0.0 {
            true => lo.to_bits() >> self.shift,
            false => 0,
        };
        let hi_key = hi.to_bits() >> self.shift;
        if hi_key < self.first_fine || lo_key > last_fine {
            return None;
        }
        let cell = |key: u32| self.cell_of_fine[(key - self.first_fine) as usize] as usize;
        Some((
            cell(lo_key.max(self.first_fine)),
            cell(hi_key.min(last_fine)),
        ))
    }

    /// Every fragment of a peptide in `candidates` falling in one of `windows`, once per
    /// window it falls in: the same matches `page_search_batch` finds in pages. Rebuilds
    /// each fragment of the window's cells exactly and compares the m/z (the reference
    /// for `matched_candidates`).
    pub fn search(
        &self,
        candidates: Range<usize>,
        windows: &[(f32, Tolerance)],
        mut on_match: impl FnMut(Theoretical),
    ) {
        let candidates = candidates.start..candidates.end.min(self.rows);
        if candidates.is_empty() {
            return;
        }
        GRID_SCRATCH.with(|scratch| {
            let searched = &mut scratch.borrow_mut().windows;
            self.cell_windows(windows, searched);
            for block in self.blocks_holding(&candidates) {
                let (groups, directory) = self.block(block);
                for slab in 0..self.slabs {
                    let starts = self.slab_starts(directory, slab);
                    for window in searched.iter() {
                        for cell in window.first..=window.last {
                            for fragment in starts[cell] as usize..starts[cell + 1] as usize {
                                let theoretical = self.fragment(block, groups, cell, fragment);
                                let mz = theoretical.fragment_mz;
                                let row = theoretical.peptide_index.0 as usize;
                                if mz >= window.lo && mz <= window.hi && candidates.contains(&row) {
                                    on_match(theoretical);
                                }
                            }
                        }
                    }
                }
            }
        });
    }

    /// Appends to `matched`, for every fragment of a peptide in `candidates` falling in
    /// one of `windows` (once per window it falls in), the peptide's offset from
    /// `candidates.start`: the matches of `search`, in another order. With an RT axis and
    /// `rt` (the predicted RTs compatible with the spectrum), rows in RT slabs wholly
    /// outside `rt` are skipped; rows in the slabs read may still lie outside it, so the
    /// caller's own RT check stays.
    pub fn matched_candidates(
        &self,
        candidates: Range<usize>,
        windows: &[(f32, Tolerance)],
        rt: Option<(f32, f32)>,
        matched: &mut Vec<u32>,
    ) {
        let visits = self.slab_visits(rt);
        let candidates = candidates.start..candidates.end.min(self.rows);
        if candidates.is_empty() {
            return;
        }
        GRID_SCRATCH.with(|scratch| {
            let searched = &mut scratch.borrow_mut().windows;
            self.cell_windows(windows, searched);
            #[cfg(target_arch = "x86_64")]
            if std::arch::is_x86_feature_detected!("avx2")
                && std::arch::is_x86_feature_detected!("popcnt")
            {
                // SAFETY: AVX2 and POPCNT support checked at run time.
                unsafe { self.matched_candidates_avx2(&candidates, searched, &visits, matched) };
                return;
            }
            self.matched_candidates_scalar(&candidates, searched, &visits, matched);
        });
    }

    fn matched_candidates_scalar(
        &self,
        candidates: &Range<usize>,
        searched: &[CellWindow],
        visits: &[Range<usize>; 2],
        matched: &mut Vec<u32>,
    ) {
        for block in self.blocks_holding(candidates) {
            let (groups, directory) = self.block(block);
            let (local, offset_of_local) = self.local_candidates(block, candidates);
            for slab in visits.iter().flat_map(Clone::clone) {
                let starts = self.slab_starts(directory, slab);
                for window in searched {
                    for (run, residual_lo, residual_hi) in window.runs(starts) {
                        for fragment in run {
                            let group = &groups[fragment / LANES];
                            let residual = group.residuals[fragment % LANES];
                            let row = group.rows[fragment % LANES] as u32;
                            if residual >= residual_lo
                                && residual <= residual_hi
                                && local.contains(&row)
                            {
                                matched.push(row.wrapping_add(offset_of_local));
                            }
                        }
                    }
                }
            }
        }
    }

    /// `matched_candidates_scalar`, 16 lanes at a time without branches: residual and
    /// row masks, lanes outside the run masked by position, then each half group's
    /// matching offsets left-packed into `matched`.
    #[cfg(target_arch = "x86_64")]
    #[target_feature(enable = "avx2,popcnt")]
    unsafe fn matched_candidates_avx2(
        &self,
        candidates: &Range<usize>,
        searched: &[CellWindow],
        visits: &[Range<usize>; 2],
        matched: &mut Vec<u32>,
    ) {
        use std::arch::x86_64::*;
        for block in self.blocks_holding(candidates) {
            let (groups, directory) = self.block(block);
            let (local, offset_of_local) = self.local_candidates(block, candidates);
            let first_local = _mm256_set1_epi16(local.start as u16 as i16);
            let last_local = _mm256_set1_epi16((local.end - local.start - 1) as u16 as i16);
            let offset_of_local = _mm256_set1_epi32(offset_of_local as i32);
            for slab in visits.iter().flat_map(Clone::clone) {
                let starts = self.slab_starts(directory, slab);
                for window in searched {
                    for (run, residual_lo, residual_hi) in window.runs(starts) {
                        if run.is_empty() {
                            continue;
                        }
                        let lower = _mm256_set1_epi16(residual_lo as i16);
                        let upper = _mm256_set1_epi16(residual_hi as i16);
                        let groups_of_run = run.start / LANES..(run.end + LANES - 1) / LANES;
                        matched.reserve(groups_of_run.len() * LANES);
                        let out = matched.as_mut_ptr();
                        let mut len = matched.len();
                        for g in groups_of_run {
                            let group = &groups[g];
                            let residuals =
                                _mm256_load_si256(group.residuals.as_ptr() as *const __m256i);
                            let rows = _mm256_load_si256(group.rows.as_ptr() as *const __m256i);
                            let in_window = _mm256_and_si256(
                                _mm256_cmpeq_epi16(_mm256_max_epu16(residuals, lower), residuals),
                                _mm256_cmpeq_epi16(_mm256_min_epu16(residuals, upper), residuals),
                            );
                            let rows_from_first = _mm256_sub_epi16(rows, first_local);
                            let is_candidate = _mm256_cmpeq_epi16(
                                _mm256_min_epu16(rows_from_first, last_local),
                                rows_from_first,
                            );
                            let both = _mm256_and_si256(in_window, is_candidate);
                            let bytes = _mm256_movemask_epi8(_mm256_packs_epi16(
                                both,
                                _mm256_setzero_si256(),
                            )) as u32;
                            let in_run = lanes_in(g, &run);
                            let low = (bytes & in_run & 0xFF) as usize;
                            let high = ((bytes >> 16) & (in_run >> HALF) & 0xFF) as usize;
                            let offsets_low = _mm256_add_epi32(
                                _mm256_cvtepu16_epi32(_mm256_castsi256_si128(rows)),
                                offset_of_local,
                            );
                            let offsets_high = _mm256_add_epi32(
                                _mm256_cvtepu16_epi32(_mm256_extracti128_si256(rows, 1)),
                                offset_of_local,
                            );
                            // SAFETY: `reserve` left room for 16 more per group of the run.
                            _mm256_storeu_si256(
                                out.add(len) as *mut __m256i,
                                _mm256_permutevar8x32_epi32(
                                    offsets_low,
                                    _mm256_loadu_si256(LEFT_PACK[low].as_ptr() as *const __m256i),
                                ),
                            );
                            len += low.count_ones() as usize;
                            _mm256_storeu_si256(
                                out.add(len) as *mut __m256i,
                                _mm256_permutevar8x32_epi32(
                                    offsets_high,
                                    _mm256_loadu_si256(LEFT_PACK[high].as_ptr() as *const __m256i),
                                ),
                            );
                            len += high.count_ones() as usize;
                        }
                        // SAFETY: the first `len` slots are written.
                        matched.set_len(len);
                    }
                }
            }
        }
    }
}

/// Upper bounds of `rt.slabs` RT slabs holding about equally many rows (from a sample of
/// the finite predictions); fewer when predictions repeat.
fn rt_slab_bounds(rt: &RtSlabbing) -> Vec<f32> {
    let step = (rt.predicted.len() / 100_000).max(1);
    let mut sample: Vec<f32> = (0..rt.predicted.len())
        .step_by(step)
        .map(|row| rt.predicted.get(row))
        .filter(|value| value.is_finite())
        .collect();
    sample.sort_unstable_by(f32::total_cmp);
    let mut bounds: Vec<f32> = (1..rt.slabs)
        .filter_map(|slab| sample.get(slab * sample.len() / rt.slabs).copied())
        .collect();
    bounds.dedup_by(|a, b| a.to_bits() == b.to_bits());
    bounds
}

/// The lanes of group `group` whose fragments lie in `run`, as a 16-bit mask.
#[inline(always)]
fn lanes_in(group: usize, run: &Range<usize>) -> u32 {
    let group_start = group * LANES;
    let first = run.start.saturating_sub(group_start).min(LANES) as u32;
    let end = (run.end - group_start).min(LANES) as u32;
    ((1u32 << end) - 1) & !((1u32 << first) - 1)
}

/// Splits fine keys `first_fine..` (with `histogram[i]` fragments at `first_fine + i`)
/// into cells of consecutive keys: a cell closes once it holds `target` fragments or
/// spans `width_cap` keys. Returns each key's cell and each cell's first key.
fn equal_count_cells(
    histogram: &[u32],
    first_fine: u32,
    target: u64,
    width_cap: usize,
) -> (Vec<u16>, Vec<u32>) {
    let mut cell_of_fine = Vec::with_capacity(histogram.len());
    let mut cell_first_fine = Vec::new();
    let mut held = 0u64;
    let mut width = 0;
    for (offset, &count) in histogram.iter().enumerate() {
        if cell_first_fine.is_empty() || held >= target || width == width_cap {
            cell_first_fine.push(first_fine + offset as u32);
            held = 0;
            width = 0;
        }
        let cell = cell_first_fine.len() - 1;
        assert!(
            cell <= u16::MAX as usize,
            "more than 65,536 grid cells: raise SAGE_GRID_CELL_FRAGMENTS"
        );
        cell_of_fine.push(cell as u16);
        held += u64::from(count);
        width += 1;
    }
    (cell_of_fine, cell_first_fine)
}

/// `values` cut at `starts` (ascending, first 0, last `values.len()`).
fn split_by_starts<'a, T>(mut values: &'a mut [T], starts: &[usize]) -> Vec<&'a mut [T]> {
    let mut pieces = Vec::with_capacity(starts.len().saturating_sub(1));
    for bounds in starts.windows(2) {
        let (piece, rest) = values.split_at_mut(bounds[1] - bounds[0]);
        pieces.push(piece);
        values = rest;
    }
    pieces
}

#[cfg(test)]
mod test {
    use super::*;

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }

        fn below(&mut self, bound: u64) -> u64 {
            self.next() % bound
        }
    }

    /// Per row 0–9 fragments over 100–2100 m/z, some repeated within and across rows.
    fn random_fragments(rows: usize, random: &mut Lcg) -> Vec<Vec<f32>> {
        let mut by_row: Vec<Vec<f32>> = Vec::new();
        for row in 0..rows {
            let mut fragments = Vec::new();
            for _ in 0..random.below(10) {
                let mz = match (row > 0 && random.below(4) == 0, by_row.last()) {
                    (true, Some(previous)) if !previous.is_empty() => {
                        previous[random.below(previous.len() as u64) as usize]
                    }
                    _ => 100.0 + random.below(2_000_000) as f32 / 1000.0,
                };
                fragments.push(mz);
                if random.below(8) == 0 {
                    fragments.push(mz);
                }
            }
            by_row.push(fragments);
        }
        by_row
    }

    fn sorted(mut offsets: Vec<u32>) -> Vec<u32> {
        offsets.sort_unstable();
        offsets
    }

    #[test]
    fn matched_candidates_equal_the_exact_search_in_every_lane_width() {
        let mut random = Lcg(7);
        let rows = 37;
        let by_row = random_fragments(rows, &mut random);
        let all: Vec<f32> = by_row.iter().flatten().copied().collect();
        let mut windows: Vec<(f32, Tolerance)> = (0..40)
            .map(|_| {
                let mz = all[random.below(all.len() as u64) as usize];
                (mz, Tolerance::Ppm(-20.0, 20.0))
            })
            .collect();
        // Windows starting, ending or both exactly on a fragment's m/z.
        for _ in 0..10 {
            let mz = all[random.below(all.len() as u64) as usize];
            windows.extend([
                (mz, Tolerance::Da(0.0, 0.0)),
                (mz, Tolerance::Da(-0.3, 0.0)),
                (mz, Tolerance::Da(0.0, 0.3)),
            ]);
        }
        windows.extend([
            (512.0, Tolerance::Da(-0.5, 0.5)),
            (1024.0, Tolerance::Da(0.0, 0.0)),
            (300.0, Tolerance::Da(-250.0, 250.0)),
            (1000.0, Tolerance::Da(-2000.0, 2000.0)),
            (50.0, Tolerance::Da(-1.0, 1.0)),
            (5000.0, Tolerance::Da(-1.0, 1.0)),
        ]);
        let candidate_ranges = [0..rows, 0..1, 5..6, 3..30, 8..16, 16..37, 36..37, 10..10];
        for rows_per_block in [1, 3, 8, 16, 64] {
            for (fine_bits, cell_fragments) in [(7, 1), (10, 4), (14, 1), (14, 16), (16, 1000)] {
                let shape = GridShape {
                    rows_per_block,
                    fine_bits,
                    cell_fragments,
                };
                let grid = FragmentGrid::build(rows, shape, None, |row, out| {
                    out.extend(by_row[row].iter().map(|&fragment_mz| Theoretical {
                        peptide_index: PeptideIx(row as u32),
                        fragment_mz,
                    }))
                });
                assert_eq!(grid.len(), all.len());
                assert_cells_partition_the_fine_keys(&grid);
                let mut stored: Vec<(u32, u32)> = grid
                    .fragments()
                    .map(|f| (f.peptide_index.0, f.fragment_mz.to_bits()))
                    .collect();
                let mut generated: Vec<(u32, u32)> = by_row
                    .iter()
                    .enumerate()
                    .flat_map(|(row, mzs)| mzs.iter().map(move |mz| (row as u32, mz.to_bits())))
                    .collect();
                stored.sort_unstable();
                generated.sort_unstable();
                assert_eq!(stored, generated, "{shape:?}: residuals rebuild every m/z");
                let mut searched = Vec::new();
                grid.cell_windows(&windows, &mut searched);
                for candidates in candidate_ranges.clone() {
                    let mut exact = Vec::new();
                    grid.search(candidates.clone(), &windows, |fragment| {
                        exact.push(fragment.peptide_index.0 - candidates.start as u32)
                    });
                    let mut dispatched = Vec::new();
                    grid.matched_candidates(candidates.clone(), &windows, None, &mut dispatched);
                    assert_eq!(
                        sorted(dispatched),
                        sorted(exact.clone()),
                        "{shape:?} {candidates:?}"
                    );
                    if candidates.is_empty() {
                        continue;
                    }
                    let mut scalar = Vec::new();
                    let all_slabs = grid.slab_visits(None);
                    grid.matched_candidates_scalar(&candidates, &searched, &all_slabs, &mut scalar);
                    assert_eq!(
                        sorted(scalar),
                        sorted(exact.clone()),
                        "{shape:?} {candidates:?}"
                    );
                    #[cfg(target_arch = "x86_64")]
                    if std::arch::is_x86_feature_detected!("avx2")
                        && std::arch::is_x86_feature_detected!("popcnt")
                    {
                        let mut vector = Vec::new();
                        // SAFETY: AVX2 and POPCNT support checked at run time.
                        unsafe {
                            grid.matched_candidates_avx2(
                                &candidates,
                                &searched,
                                &all_slabs,
                                &mut vector,
                            )
                        };
                        assert_eq!(sorted(vector), sorted(exact), "{shape:?} {candidates:?}");
                    }
                }
            }
        }
    }

    /// Cells are runs of consecutive fine keys, numbered in order, each within the width
    /// cap, and every fragment sits in its own fine key's cell.
    fn assert_cells_partition_the_fine_keys(grid: &FragmentGrid) {
        let cap = 1 << (23 - grid.shift - MIN_FINE_BITS);
        assert_eq!(grid.cell_of_fine.first().copied().unwrap_or(0), 0);
        for (offset, pair) in grid.cell_of_fine.windows(2).enumerate() {
            assert!(
                pair[1] == pair[0] || pair[1] == pair[0] + 1,
                "fine key {offset}"
            );
        }
        assert_eq!(grid.cell_first_fine.len(), grid.cells);
        for (cell, &first) in grid.cell_first_fine.iter().enumerate() {
            let offset = (first - grid.first_fine) as usize;
            assert_eq!(grid.cell_of_fine[offset] as usize, cell);
            let width = grid.cell_of_fine[offset..]
                .iter()
                .take_while(|&&c| c as usize == cell)
                .count();
            assert!(width <= cap, "cell {cell} spans {width} fine keys");
        }
        for block in 0..grid.block_groups.len() - 1 {
            let (groups, starts) = grid.block(block);
            for cell in 0..grid.cells {
                for fragment in starts[cell] as usize..starts[cell + 1] as usize {
                    let mz = grid.fragment(block, groups, cell, fragment).fragment_mz;
                    let key = (mz.to_bits() >> grid.shift) - grid.first_fine;
                    assert_eq!(grid.cell_of_fine[key as usize] as usize, cell);
                }
            }
        }
    }

    /// `small` minus `large`, as multisets: empty when `small` is contained in `large`.
    /// A census over every `every`-th block (and the last) finds exactly what counting
    /// every block finds; with a fragment outside the sampled range in an unsampled
    /// block it falls back to the exact census.
    #[test]
    fn sampled_census_finds_what_the_exact_census_finds() {
        let mut random = Lcg(41);
        let rows = 200;
        let mut by_row = random_fragments(rows, &mut random);
        let all: Vec<f32> = by_row.iter().flatten().copied().collect();
        let windows: Vec<(f32, Tolerance)> = (0..40)
            .map(|_| {
                let mz = all[random.below(all.len() as u64) as usize];
                (mz, Tolerance::Ppm(-30.0, 30.0))
            })
            .chain([(1000.0, Tolerance::Da(-2000.0, 2000.0))])
            .collect();
        let grids = |by_row: &Vec<Vec<f32>>, shape: GridShape, every: usize| {
            let fragments_of = |row: usize, out: &mut Vec<Theoretical>| {
                out.extend(by_row[row].iter().map(|&fragment_mz| Theoretical {
                    peptide_index: PeptideIx(row as u32),
                    fragment_mz,
                }))
            };
            let fragment_count = |row: usize| by_row[row].len();
            let sampled = FragmentGrid::build_sampled(
                rows,
                shape,
                None,
                &fragments_of,
                Some(&fragment_count),
                every,
            );
            (sampled, FragmentGrid::build(rows, shape, None, fragments_of))
        };
        let fragment_set = |grid: &FragmentGrid| {
            let mut fragments: Vec<(u32, u32)> = grid
                .fragments()
                .map(|f| (f.peptide_index.0, f.fragment_mz.to_bits()))
                .collect();
            fragments.sort_unstable();
            fragments
        };
        for rows_per_block in [1, 4, 16] {
            for (fine_bits, cell_fragments) in [(10, 1), (14, 4)] {
                let shape = GridShape {
                    rows_per_block,
                    fine_bits,
                    cell_fragments,
                };
                for every in [2, 3, 5] {
                    let (sampled, exact) = grids(&by_row, shape, every);
                    assert_cells_partition_the_fine_keys(&sampled);
                    assert_eq!(fragment_set(&sampled), fragment_set(&exact));
                    for candidates in [0..rows, 3..150, 77..78] {
                        let (mut a, mut b) = (Vec::new(), Vec::new());
                        sampled.matched_candidates(candidates.clone(), &windows, None, &mut a);
                        exact.matched_candidates(candidates.clone(), &windows, None, &mut b);
                        assert_eq!(sorted(a), sorted(b), "{shape:?}, every {every}");
                    }
                }
            }
        }
        // Extremes in unsampled blocks (rows 5 and 9, blocks 1 and 2 at R 4, every 4).
        by_row[5].push(30_000.0);
        by_row[9].push(17.0);
        let shape = GridShape {
            rows_per_block: 4,
            fine_bits: 10,
            cell_fragments: 1,
        };
        let (sampled, exact) = grids(&by_row, shape, 4);
        assert_eq!(sampled.cell_first_fine, exact.cell_first_fine, "fell back to the exact census");
        assert_eq!(fragment_set(&sampled), fragment_set(&exact));
    }

    fn missing_from(small: &[u32], large: &[u32]) -> Vec<u32> {
        let mut large = sorted(large.to_vec());
        let mut missing = Vec::new();
        for value in sorted(small.to_vec()) {
            match large.binary_search(&value) {
                Ok(at) => {
                    large.remove(at);
                }
                Err(_) => missing.push(value),
            }
        }
        missing
    }

    #[test]
    fn rt_slabs_keep_every_match_sage_keeps_and_add_none() {
        let mut random = Lcg(29);
        let rows = 61;
        let by_row = random_fragments(rows, &mut random);
        let predicted: Vec<f64> = (0..rows)
            .map(|_| match random.below(50) {
                0..=5 => f64::NAN,
                6 => f64::INFINITY,
                _ => random.below(10_000) as f64 / 100.0,
            })
            .collect();
        let all: Vec<f32> = by_row.iter().flatten().copied().collect();
        let windows: Vec<(f32, Tolerance)> = (0..40)
            .map(|_| {
                (
                    all[random.below(all.len() as u64) as usize],
                    Tolerance::Ppm(-30.0, 30.0),
                )
            })
            .chain([(1000.0, Tolerance::Da(-2000.0, 2000.0))])
            .collect();
        let rt_windows = [
            None,
            Some((10.0, 20.0)),
            Some((0.0, 100.0)),
            Some((55.5, 55.6)),
            Some((-5.0, 3.0)),
            Some((99.0, 150.0)),
            Some((30.0, 20.0)),
            Some((f32::NAN, 10.0)),
        ];
        for rows_per_block in [1, 7, 64] {
            for slabs in [1, 2, 5, 13] {
                let shape = GridShape {
                    rows_per_block,
                    fine_bits: 10,
                    cell_fragments: 2,
                };
                let rt = RtSlabbing {
                    predicted: PredictedValues::new(&predicted),
                    slabs,
                };
                let grid = FragmentGrid::build(rows, shape, Some(rt), |row, out| {
                    out.extend(by_row[row].iter().map(|&fragment_mz| Theoretical {
                        peptide_index: PeptideIx(row as u32),
                        fragment_mz,
                    }))
                });
                assert_eq!(grid.len(), all.len());
                for candidates in [0..rows, 3..40, 20..21] {
                    let mut unpruned = Vec::new();
                    grid.search(candidates.clone(), &windows, |fragment| {
                        unpruned.push(fragment.peptide_index.0 - candidates.start as u32)
                    });
                    // Windows starting or ending exactly on a slab bound, where rows can sit.
                    let on_bounds = grid.rt_bounds.iter().flat_map(|&bound| {
                        [
                            Some((bound, bound)),
                            Some((bound, bound + 1.0)),
                            Some((bound - 1.0, bound)),
                        ]
                    });
                    for rt_window in rt_windows.into_iter().chain(on_bounds) {
                        let kept: Vec<u32> = unpruned
                            .iter()
                            .copied()
                            .filter(|&offset| {
                                let rt = predicted[candidates.start + offset as usize] as f32;
                                match rt_window {
                                    None => true,
                                    Some((lo, hi)) => rt >= lo && rt <= hi,
                                }
                            })
                            .collect();
                        let mut pruned = Vec::new();
                        grid.matched_candidates(
                            candidates.clone(),
                            &windows,
                            rt_window,
                            &mut pruned,
                        );
                        let label = format!(
                            "R {rows_per_block}, slabs {slabs}, {candidates:?}, {rt_window:?}"
                        );
                        assert!(
                            missing_from(&pruned, &unpruned).is_empty(),
                            "{label}: extra"
                        );
                        assert!(missing_from(&kept, &pruned).is_empty(), "{label}: lost");
                        if rt_window.is_none() {
                            assert_eq!(sorted(pruned.clone()), sorted(unpruned.clone()), "{label}");
                        }
                        let mut searched = Vec::new();
                        grid.cell_windows(&windows, &mut searched);
                        let visits = grid.slab_visits(rt_window);
                        let mut scalar = Vec::new();
                        grid.matched_candidates_scalar(
                            &candidates,
                            &searched,
                            &visits,
                            &mut scalar,
                        );
                        assert_eq!(sorted(scalar), sorted(pruned.clone()), "{label}: scalar");
                    }
                }
            }
        }
    }

    #[test]
    fn equal_count_cells_close_on_count_or_width() {
        let (cell_of_fine, first) = equal_count_cells(&[5, 0, 0, 7, 1, 1, 1, 1, 0, 9], 100, 6, 3);
        assert_eq!(cell_of_fine, vec![0, 0, 0, 1, 2, 2, 2, 3, 3, 3]);
        assert_eq!(first, vec![100, 103, 104, 107]);
    }

    #[test]
    fn left_pack_lists_set_lanes_in_order() {
        for (mask, packed) in LEFT_PACK.iter().enumerate() {
            let lanes: Vec<u32> = (0..LANES as u32)
                .filter(|lane| mask & (1 << lane) != 0)
                .collect();
            assert_eq!(&packed[..lanes.len()], &lanes[..]);
        }
    }

    #[test]
    fn lanes_in_masks_the_run_by_position() {
        assert_eq!(lanes_in(0, &(0..16)), 0xFFFF);
        assert_eq!(lanes_in(0, &(3..16)), 0xFFF8);
        assert_eq!(lanes_in(1, &(3..20)), 0x000F);
        assert_eq!(lanes_in(1, &(17..18)), 0x0002);
        assert_eq!(lanes_in(2, &(0..100)), 0xFFFF);
    }

    #[test]
    fn fragment_group_is_one_cache_line() {
        assert_eq!(std::mem::size_of::<FragmentGroup>(), 64);
        assert_eq!(std::mem::align_of::<FragmentGroup>(), 64);
    }
}
