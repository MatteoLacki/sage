//! What fragment generation and scoring read of a peptide, so they run the same on
//! a sage [`Peptide`] or straight on a [`CompactPeptides`](crate::compact::CompactPeptides)
//! row ([`CompactView`](crate::compact::CompactView)), without building a `Peptide`.

use crate::ion_series::{IonSeries, Kind};
use crate::peptide::Peptide;

/// A peptide's residues, their modification masses, terminal modifications, mass
/// and the attributes features report.
pub trait PeptideView {
    /// Number of residues.
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    fn residue(&self, index: usize) -> u8;

    /// The (static or variable) modification mass on residue `index`, 0 when none.
    fn modification(&self, index: usize) -> f32;

    fn nterm(&self) -> Option<f32>;

    fn cterm(&self) -> Option<f32>;

    /// Monoisotopic mass, modifications included.
    fn monoisotopic(&self) -> f32;

    fn decoy(&self) -> bool;

    fn missed_cleavages(&self) -> u8;

    /// 1 for targets, -1 for decoys.
    fn label(&self) -> i32 {
        match self.decoy() {
            true => -1,
            false => 1,
        }
    }
}

impl PeptideView for Peptide {
    fn len(&self) -> usize {
        self.sequence.len()
    }

    fn residue(&self, index: usize) -> u8 {
        self.sequence[index]
    }

    fn modification(&self, index: usize) -> f32 {
        self.modifications[index]
    }

    fn nterm(&self) -> Option<f32> {
        self.nterm
    }

    fn cterm(&self) -> Option<f32> {
        self.cterm
    }

    fn monoisotopic(&self) -> f32 {
        self.monoisotopic
    }

    fn decoy(&self) -> bool {
        self.decoy
    }

    fn missed_cleavages(&self) -> u8 {
        self.missed_cleavages
    }
}

/// The first way `a` and `b` differ, as fragment generation and scoring see them:
/// residues, modification bits, terminal modifications, mass, decoy, missed
/// cleavages and every ion of every kind, bit for bit. `None` when they are the same.
pub fn first_difference(a: &impl PeptideView, b: &impl PeptideView) -> Option<String> {
    if a.len() != b.len() {
        return Some(format!("length {} vs {}", a.len(), b.len()));
    }
    for index in 0..a.len() {
        if a.residue(index) != b.residue(index) {
            return Some(format!("residue {index}"));
        }
        if a.modification(index).to_bits() != b.modification(index).to_bits() {
            return Some(format!(
                "modification {index}: {} vs {}",
                a.modification(index),
                b.modification(index)
            ));
        }
    }
    let terminal_bits = |mass: Option<f32>| mass.map(f32::to_bits);
    if terminal_bits(a.nterm()) != terminal_bits(b.nterm()) {
        return Some(format!("nterm {:?} vs {:?}", a.nterm(), b.nterm()));
    }
    if terminal_bits(a.cterm()) != terminal_bits(b.cterm()) {
        return Some(format!("cterm {:?} vs {:?}", a.cterm(), b.cterm()));
    }
    if a.monoisotopic().to_bits() != b.monoisotopic().to_bits() {
        return Some(format!("mass {} vs {}", a.monoisotopic(), b.monoisotopic()));
    }
    if a.decoy() != b.decoy() || a.missed_cleavages() != b.missed_cleavages() {
        return Some("decoy or missed cleavages".into());
    }
    for kind in [Kind::A, Kind::B, Kind::C, Kind::X, Kind::Y, Kind::Z] {
        if ion_bits(a, kind) != ion_bits(b, kind) {
            return Some(format!("{kind:?} ions"));
        }
    }
    None
}

fn ion_bits(peptide: &impl PeptideView, kind: Kind) -> Vec<u32> {
    IonSeries::new(peptide, kind)
        .map(|ion| ion.monoisotopic_mass.to_bits())
        .collect()
}
