//! Cutting mass-sorted peptide rows into slices whose fragment indices are built
//! and searched one at a time. A query searched in its slice finds exactly the rows
//! it would find in the whole table, including the padding row `binary_search_slice`
//! takes below a window, so scoring cannot tell the two apart.

use std::ops::Range;

pub struct Slice {
    /// Rows this slice's index holds: the core, one padding row below it, and every
    /// row lighter than the next core's first mass plus the widest query window.
    pub rows: Range<usize>,
    /// The rows whose masses decide which queries belong to this slice.
    pub core: Range<usize>,
}

pub struct SlicePlan {
    pub slices: Vec<Slice>,
    /// Mass of each core's first row, strictly increasing.
    core_start_masses: Vec<f32>,
}

impl SlicePlan {
    /// Cut rows `0..rows` (ascending by `mass(row)`) into cores costing about `budget`
    /// each, splitting only between rows of different masses, and extend each core so a
    /// query window up to `widest_window` wide that starts inside it fits whole.
    pub fn new(
        rows: usize,
        mass: impl Fn(usize) -> f32,
        cost: impl Fn(usize) -> usize,
        budget: usize,
        widest_window: f32,
    ) -> Self {
        let mut cores = Vec::new();
        let mut core_start = 0;
        let mut core_cost = 0;
        let mut previous_mass = f32::NEG_INFINITY;
        for row in 0..rows {
            let row_mass = mass(row);
            let starts_new_mass = row > core_start && previous_mass < row_mass;
            if core_cost >= budget && starts_new_mass {
                cores.push(core_start..row);
                core_start = row;
                core_cost = 0;
            }
            core_cost += cost(row);
            previous_mass = row_mass;
        }
        cores.push(core_start..rows);

        let slices = cores
            .iter()
            .map(|core| {
                let end = match core.end < rows {
                    false => rows,
                    true => {
                        let limit = mass(core.end) + widest_window;
                        partition_point(rows, |row| mass(row) < limit)
                    }
                };
                Slice {
                    rows: core.start.saturating_sub(1)..end,
                    core: core.clone(),
                }
            })
            .collect();
        let core_start_masses = cores
            .iter()
            .filter(|core| core.start < rows)
            .map(|core| mass(core.start))
            .collect();
        Self {
            slices,
            core_start_masses,
        }
    }

    /// The slice for a query whose lowest candidate mass is `window_low`: the one
    /// whose core holds that mass (the first slice for anything lighter).
    pub fn slice_for(&self, window_low: f32) -> usize {
        self.core_start_masses
            .partition_point(|&start| start <= window_low)
            .saturating_sub(1)
    }
}

/// The rows a slice must hold for queries whose windows reach from `low` (their lowest
/// window start) to `high` (their highest window end): any of those windows finds the
/// same rows in them as in all rows `0..rows` (ascending by `mass(row)`), the padding
/// row below it included.
pub fn rows_reaching(
    rows: usize,
    mass: impl Fn(usize) -> f32,
    low: f32,
    high: f32,
) -> Range<usize> {
    let start = partition_point(rows, |row| mass(row) < low).saturating_sub(1);
    let end = partition_point(rows, |row| mass(row) <= high);
    start..end.max(start)
}

/// `slice::partition_point` over rows `0..rows`: the first row for which `before` is
/// false, `before` being true for a prefix of the rows.
fn partition_point(rows: usize, before: impl Fn(usize) -> bool) -> usize {
    let (mut low, mut high) = (0, rows);
    while low < high {
        let middle = low + (high - low) / 2;
        match before(middle) {
            true => low = middle + 1,
            false => high = middle,
        }
    }
    low
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::binary_search_slice;

    /// Deterministic pseudo-random numbers in [0, 1).
    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> f32 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 40) as f32 / (1u64 << 24) as f32
        }
    }

    fn sorted_masses_with_duplicates(random: &mut Lcg, count: usize) -> Vec<f32> {
        let mut masses: Vec<f32> = (0..count)
            .map(|_| 600.0 + (random.next() * 400.0 * 64.0).floor() / 64.0)
            .collect();
        masses.sort_by(f32::total_cmp);
        masses
    }

    fn global_range(masses: &[f32], low: f32, high: f32) -> (usize, usize) {
        binary_search_slice(masses, |m, bound| m.total_cmp(bound), low, high)
    }

    #[test]
    fn every_window_finds_its_global_rows_inside_its_slice() {
        let mut random = Lcg(7);
        for (count, budget) in [(2_000, 37), (5_000, 400), (300, 1), (1_000, 10_000)] {
            let masses = sorted_masses_with_duplicates(&mut random, count);
            let widest_window = 3.1;
            let plan = SlicePlan::new(
                masses.len(),
                |row| masses[row],
                |_| 1,
                budget,
                widest_window,
            );

            let mut windows: Vec<(f32, f32)> = (0..3_000)
                .map(|_| {
                    let low = 590.0 + random.next() * 420.0;
                    (low, low + random.next() * widest_window)
                })
                .collect();
            for core in plan.slices.iter().map(|slice| &slice.core) {
                if let Some(&start) = masses.get(core.start) {
                    windows.push((start, start + widest_window));
                    windows.push((start, start));
                }
            }

            for (low, high) in windows {
                let slice = &plan.slices[plan.slice_for(low)];
                let (left, right) = global_range(&masses, low, high);
                let rows = &masses[slice.rows.clone()];
                let (local_left, local_right) = global_range(rows, low, high);
                assert_eq!(
                    (
                        local_left + slice.rows.start,
                        local_right + slice.rows.start
                    ),
                    (left, right),
                    "window {low}..{high} in slice {:?}",
                    slice.rows
                );
            }
        }
    }

    #[test]
    fn windows_find_their_global_rows_inside_the_rows_their_slice_reaches() {
        let mut random = Lcg(13);
        for (count, budget) in [(2_000, 37), (5_000, 400), (1_000, 10_000)] {
            let masses = sorted_masses_with_duplicates(&mut random, count);
            let widest_window = 3.1;
            let plan = SlicePlan::new(
                masses.len(),
                |row| masses[row],
                |_| 1,
                budget,
                widest_window,
            );
            let mut windows_by_slice = vec![Vec::new(); plan.slices.len()];
            for _ in 0..3_000 {
                let low = 590.0 + random.next() * 420.0;
                let window = (low, low + random.next() * widest_window);
                windows_by_slice[plan.slice_for(low)].push(window);
            }
            for (slice, windows) in plan.slices.iter().zip(&windows_by_slice) {
                if windows.is_empty() {
                    continue;
                }
                let low = windows.iter().map(|w| w.0).fold(f32::INFINITY, f32::min);
                let high = windows
                    .iter()
                    .map(|w| w.1)
                    .fold(f32::NEG_INFINITY, f32::max);
                let reach = rows_reaching(masses.len(), |row| masses[row], low, high);
                assert!(slice.rows.start <= reach.start && reach.end <= slice.rows.end);
                for &(low, high) in windows {
                    let (left, right) = global_range(&masses, low, high);
                    let (local_left, local_right) = global_range(&masses[reach.clone()], low, high);
                    assert_eq!(
                        (local_left + reach.start, local_right + reach.start),
                        (left, right)
                    );
                }
            }
        }
    }

    #[test]
    fn cores_partition_the_rows_at_mass_changes() {
        let mut random = Lcg(11);
        let masses = sorted_masses_with_duplicates(&mut random, 4_000);
        let plan = SlicePlan::new(masses.len(), |row| masses[row], |_| 1, 50, 3.1);
        let mut next = 0;
        for slice in &plan.slices {
            assert_eq!(slice.core.start, next);
            assert!(slice.core.start < slice.core.end);
            if slice.core.start > 0 {
                assert!(masses[slice.core.start - 1] < masses[slice.core.start]);
            }
            next = slice.core.end;
        }
        assert_eq!(next, masses.len());
        assert!(plan.slices.len() > 10);
    }

    #[test]
    fn empty_rows_make_one_empty_slice() {
        let plan = SlicePlan::new(0, |_| unreachable!(), |_| 1, 10, 3.1);
        assert_eq!(plan.slices.len(), 1);
        assert_eq!(plan.slices[0].rows, 0..0);
        assert_eq!(plan.slice_for(1000.0), 0);
    }
}
