//! Per-peptide-row predictions read in place, in the layout their producer
//! (`git/featureprediction`, `predict_prosit`) stores them: one value per peptide row,
//! or one per (peptide row, precursor charge) slot of a SPECLIB, numbered as its
//! `slot_order` says ([`SlotOrder`]).

use std::ops::Range;

/// Predicted RT values, stored as `f64` and compared as `f32`.
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

/// How a SPECLIB numbers its (peptide row, precursor charge) slots: its manifest's
/// `slot_order`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlotOrder {
    /// `charge_major_v1`: slot = (charge − min_charge) · rows + row.
    ChargeMajor,
    /// `row_major_v1`: slot = row · charge_span + (charge − min_charge).
    RowMajor,
}

impl SlotOrder {
    pub fn from_token(token: &str) -> Option<Self> {
        match token {
            "charge_major_v1" => Some(Self::ChargeMajor),
            "row_major_v1" => Some(Self::RowMajor),
            _ => None,
        }
    }
}

/// The slots of a table of `rows` peptide rows at precursor charges
/// `min_charge..=max_charge`, seen from a scorer whose peptide `i` is table row
/// `first_row + i`.
#[derive(Clone, Copy, Debug)]
pub struct Slots {
    pub order: SlotOrder,
    pub rows: usize,
    pub first_row: usize,
    pub min_charge: u8,
    pub max_charge: u8,
}

impl Slots {
    pub fn charge_span(&self) -> usize {
        (self.max_charge - self.min_charge + 1) as usize
    }

    /// How many slots the whole table has.
    pub fn count(&self) -> usize {
        self.rows * self.charge_span()
    }

    /// The slot of the scorer's peptide `peptide` at `charge`; `None` for a charge
    /// outside the range.
    pub fn get(&self, peptide: usize, charge: u8) -> Option<usize> {
        if charge < self.min_charge || charge > self.max_charge {
            return None;
        }
        let row = self.first_row + peptide;
        let charge_offset = (charge - self.min_charge) as usize;
        Some(match self.order {
            SlotOrder::ChargeMajor => charge_offset * self.rows + row,
            SlotOrder::RowMajor => row * self.charge_span() + charge_offset,
        })
    }

    /// The same table, seen from a scorer whose peptide `i` is table row `first_row + i`.
    pub fn from_row(self, first_row: usize) -> Self {
        Self { first_row, ..self }
    }
}

/// One `f32` per slot (predicted IIM); NaN where there is no prediction.
#[derive(Clone, Copy, Debug)]
pub struct SlotValues<'a> {
    values: &'a [f32],
    slots: Slots,
}

impl<'a> SlotValues<'a> {
    pub fn new(values: &'a [f32], slots: Slots) -> Self {
        assert_eq!(values.len(), slots.count(), "one value per slot");
        Self { values, slots }
    }

    pub fn get(&self, peptide: usize, charge: u8) -> Option<f32> {
        self.slots.get(peptide, charge).map(|slot| self.values[slot])
    }
}

/// Fragment ranges per slot: slot `i` is `starts[i]..starts[i + 1]` of the SPECLIB's
/// arrays; an empty range is a slot without a prediction.
#[derive(Clone, Copy, Debug)]
pub struct SlotRanges<'a> {
    starts: &'a [u64],
    slots: Slots,
}

impl<'a> SlotRanges<'a> {
    pub fn new(starts: &'a [u64], slots: Slots) -> Self {
        assert_eq!(starts.len(), slots.count() + 1, "one start per slot, plus the end");
        Self { starts, slots }
    }

    /// The range of the scorer's peptide `peptide` at `charge`, if it has fragments.
    pub fn get(&self, peptide: usize, charge: u8) -> Option<(usize, usize)> {
        let slot = self.slots.get(peptide, charge)?;
        let (start, end) = (self.starts[slot] as usize, self.starts[slot + 1] as usize);
        (end > start).then_some((start, end))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn slots(order: SlotOrder) -> Slots {
        Slots {
            order,
            rows: 6,
            first_row: 0,
            min_charge: 2,
            max_charge: 4,
        }
    }

    #[test]
    fn row_major_slots() {
        let slots = slots(SlotOrder::RowMajor);
        assert_eq!(slots.get(0, 2), Some(0));
        assert_eq!(slots.get(0, 4), Some(2));
        assert_eq!(slots.get(1, 2), Some(3));
        assert_eq!(slots.get(5, 4), Some(17));
        assert_eq!(slots.from_row(4).get(1, 4), Some(17));
    }

    #[test]
    fn charge_major_slots() {
        let slots = slots(SlotOrder::ChargeMajor);
        assert_eq!(slots.get(0, 2), Some(0));
        assert_eq!(slots.get(5, 2), Some(5));
        assert_eq!(slots.get(0, 3), Some(6));
        assert_eq!(slots.get(5, 4), Some(17));
        assert_eq!(slots.from_row(4).get(1, 3), Some(11));
    }

    #[test]
    fn charges_outside_the_range_have_no_slot() {
        for order in [SlotOrder::RowMajor, SlotOrder::ChargeMajor] {
            assert_eq!(slots(order).get(0, 1), None);
            assert_eq!(slots(order).get(0, 5), None);
        }
    }

    #[test]
    fn tokens() {
        assert_eq!(SlotOrder::from_token("charge_major_v1"), Some(SlotOrder::ChargeMajor));
        assert_eq!(SlotOrder::from_token("row_major_v1"), Some(SlotOrder::RowMajor));
        assert_eq!(SlotOrder::from_token("peptide_row * charge_span + charge - min_charge"), None);
    }

    #[test]
    fn empty_ranges_are_no_prediction() {
        let slots = Slots {
            order: SlotOrder::ChargeMajor,
            rows: 2,
            first_row: 0,
            min_charge: 2,
            max_charge: 3,
        };
        // slots: (row 0, z2) 0..3, (row 1, z2) empty, (row 0, z3) 3..4, (row 1, z3) 4..6
        let starts = [0u64, 3, 3, 4, 6];
        let ranges = SlotRanges::new(&starts, slots);
        assert_eq!(ranges.get(0, 2), Some((0, 3)));
        assert_eq!(ranges.get(1, 2), None);
        assert_eq!(ranges.get(0, 3), Some((3, 4)));
        assert_eq!(ranges.get(1, 3), Some((4, 6)));
        let from_row_1 = SlotRanges::new(&starts, slots.from_row(1));
        assert_eq!(from_row_1.get(0, 3), Some((4, 6)));
    }
}
