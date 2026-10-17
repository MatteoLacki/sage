//! Per-peptide-row predictions read in place, in the layout their producer
//! (`git/featureprediction`) stores them: one value per peptide row, or per peptide
//! row and precursor charge slot (`Scorer::iim_dense_slot`).

use std::ops::Range;

/// Predicted RT or IIM values, stored as `f64` and compared as `f32`.
#[derive(Clone, Copy, Debug)]
pub struct PredictedValues<'a>(&'a [f64]);

impl<'a> PredictedValues<'a> {
    pub fn new(values: &'a [f64]) -> Self {
        Self(values)
    }

    pub fn get(&self, index: usize) -> f32 {
        self.0[index] as f32
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn slice(&self, range: Range<usize>) -> Self {
        Self(&self.0[range])
    }
}

/// Half-open ranges into the fragment-intensity cache's arrays, per (row, charge)
/// slot; a negative start marks a slot the cache has no entry for.
#[derive(Clone, Copy, Debug)]
pub struct PredictedRanges<'a> {
    start: &'a [i64],
    end: &'a [i64],
}

impl<'a> PredictedRanges<'a> {
    pub fn new(start: &'a [i64], end: &'a [i64]) -> Self {
        assert_eq!(start.len(), end.len(), "one end per start");
        Self { start, end }
    }

    /// Slot `slot`'s range, if it has one.
    pub fn get(&self, slot: usize) -> Option<(usize, usize)> {
        let start = *self.start.get(slot)?;
        if start < 0 {
            return None;
        }
        Some((start as usize, self.end[slot] as usize))
    }

    pub fn len(&self) -> usize {
        self.start.len()
    }

    pub fn is_empty(&self) -> bool {
        self.start.is_empty()
    }

    pub fn slice(&self, range: Range<usize>) -> Self {
        Self {
            start: &self.start[range.clone()],
            end: &self.end[range],
        }
    }
}
