use std::{iter::Sum, ops::Mul};

use serde::{Deserialize, Serialize};

pub const H2O: f32 = 18.010565;
pub const PROTON: f32 = 1.0072764;
pub const NEUTRON: f32 = 1.00335;
pub const NH3: f32 = 17.026548;

/// A window `[lo, hi]` on `observed − theoretical`, where "theoretical"
/// is whatever the observation is compared against: a peptide or fragment
/// mass, or an externally predicted RT/IIM. Ppm/Pct windows are relative
/// to whichever of the two values is known.
#[derive(Copy, Clone, Serialize, Deserialize, Debug, PartialEq, PartialOrd)]
#[serde(rename_all = "lowercase")]
pub enum Tolerance {
    Ppm(f32, f32),
    Pct(f32, f32),
    Da(f32, f32),
}

impl Tolerance {
    /// `(lo, hi)` as absolute deltas at `value`.
    fn absolute_deltas(&self, value: f32) -> (f32, f32) {
        match self {
            Tolerance::Ppm(lo, hi) => (value * lo / 1_000_000.0, value * hi / 1_000_000.0),
            Tolerance::Pct(lo, hi) => (value * lo / 100.0, value * hi / 100.0),
            Tolerance::Da(lo, hi) => (*lo, *hi),
        }
    }

    /// Theoretical values compatible with `observed`: `[observed − hi, observed − lo]`.
    pub fn theoretical_range(&self, observed: f32) -> (f32, f32) {
        let (lo, hi) = self.absolute_deltas(observed);
        (observed - hi, observed - lo)
    }

    /// Observed values compatible with `theoretical`: `[theoretical + lo, theoretical + hi]`.
    pub fn observed_range(&self, theoretical: f32) -> (f32, f32) {
        let (lo, hi) = self.absolute_deltas(theoretical);
        (theoretical + lo, theoretical + hi)
    }

    pub fn ppm_to_delta_mass(center: f32, ppm: f32) -> f32 {
        ppm * center / 1_000_000.0
    }
}

impl Mul<f32> for Tolerance {
    type Output = Tolerance;

    fn mul(self, rhs: f32) -> Self::Output {
        match self {
            Tolerance::Ppm(lo, hi) => Tolerance::Ppm(lo * rhs, hi * rhs),
            Tolerance::Pct(lo, hi) => Tolerance::Pct(lo * rhs, hi * rhs),
            Tolerance::Da(lo, hi) => Tolerance::Da(lo * rhs, hi * rhs),
        }
    }
}

pub const VALID_AA: [u8; 22] = [
    b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S',
    b'T', b'V', b'W', b'Y', b'U', b'O',
];

pub const MONOISOTOPIC_MASSES: [f32; 26] = [
    71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891, 113.08406, 0.0,
    128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276, 128.05858, 156.1011, 87.03203,
    101.04768, 150.95363, 99.06841, 186.07932, 0.0, 163.06332, 0.0,
];

pub const fn monoisotopic(aa: u8) -> f32 {
    if aa.is_ascii_uppercase() {
        MONOISOTOPIC_MASSES[(aa - b'A') as usize]
    } else {
        0.0
    }
}

pub const fn composition(aa: u8) -> Composition {
    match aa {
        b'A' => Composition::new(3, 2, 0),
        b'R' => Composition::new(6, 2, 0),
        b'N' => Composition::new(4, 3, 0),
        b'D' => Composition::new(4, 4, 0),
        b'C' => Composition::new(3, 2, 1),
        b'E' => Composition::new(5, 4, 0),
        b'Q' => Composition::new(5, 3, 0),
        b'G' => Composition::new(2, 2, 0),
        b'H' => Composition::new(6, 2, 0),
        b'I' => Composition::new(6, 2, 0),
        b'L' => Composition::new(6, 2, 0),
        b'K' => Composition::new(6, 2, 0),
        b'M' => Composition::new(5, 2, 1),
        b'F' => Composition::new(9, 2, 0),
        b'P' => Composition::new(5, 2, 0),
        b'S' => Composition::new(3, 3, 0),
        b'T' => Composition::new(4, 3, 0),
        b'W' => Composition::new(11, 2, 0),
        b'Y' => Composition::new(9, 3, 0),
        b'V' => Composition::new(5, 2, 0),
        b'U' => Composition::new(3, 2, 0),
        b'O' => Composition::new(12, 3, 0),
        _ => Composition::new(0, 0, 0),
    }
}

#[derive(Clone, Debug)]
pub struct Composition {
    pub carbon: u16,
    pub sulfur: u16,
}

impl Composition {
    pub const fn new(carbon: u16, _oxygen: u8, sulfur: u16) -> Self {
        Self { carbon, sulfur }
    }
}

impl Sum for Composition {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut comp = Composition::new(0, 0, 0);
        for i in iter {
            comp.carbon += i.carbon;
            // comp.oxygen += i.oxygen;
            comp.sulfur += i.sulfur;
        }
        comp
    }
}

#[cfg(test)]
mod test {
    use crate::mass::monoisotopic;

    use super::{Tolerance, VALID_AA};

    #[test]
    fn smoke() {
        for ch in VALID_AA {
            assert!(monoisotopic(ch) > 0.0);
        }
    }

    #[test]
    fn tolerances() {
        assert_eq!(
            Tolerance::Ppm(-10.0, 20.0).observed_range(1000.0),
            (999.99, 1000.02)
        );
        assert_eq!(
            Tolerance::Ppm(-10.0, 10.0).observed_range(487.0),
            (486.99513, 487.00487)
        );
        assert_eq!(
            Tolerance::Ppm(-50.0, 50.0).observed_range(1000.0),
            (999.95, 1000.05)
        );
    }

    #[test]
    fn asymmetric_windows_are_on_observed_minus_theoretical() {
        // observed − theoretical ∈ [-1, 3]
        assert_eq!(
            Tolerance::Da(-1.0, 3.0).theoretical_range(100.0),
            (97.0, 101.0)
        );
        assert_eq!(
            Tolerance::Da(-1.0, 3.0).observed_range(100.0),
            (99.0, 103.0)
        );
        assert_eq!(
            Tolerance::Ppm(-10.0, 20.0).theoretical_range(1000.0),
            (999.98, 1000.01)
        );
        assert_eq!(
            Tolerance::Pct(-1.0, 2.0).theoretical_range(100.0),
            (98.0, 101.0)
        );
        assert_eq!(
            Tolerance::Pct(-1.0, 2.0).observed_range(100.0),
            (99.0, 102.0)
        );
    }

    /// `t ∈ theoretical_range(o)` ⇔ `o ∈ observed_range(t)`: both methods
    /// describe the same window, up to the ppm²-order difference between
    /// scaling by `o` and by `t` (Da windows: exact up to f32 rounding).
    #[test]
    fn theoretical_and_observed_ranges_describe_one_window() {
        for tol in [Tolerance::Da(-0.02, 0.05), Tolerance::Ppm(-7.0, 13.0)] {
            let observed = 1000.0f32;
            let (t_lo, t_hi) = tol.theoretical_range(observed);
            for step in -30..=30 {
                let theoretical = observed + step as f32 * 0.0025;
                let (o_lo, o_hi) = tol.observed_range(theoretical);
                let margin = 1e-4;
                if theoretical > t_lo + margin && theoretical < t_hi - margin {
                    assert!(
                        observed >= o_lo && observed <= o_hi,
                        "{tol:?} t={theoretical}"
                    );
                }
                if theoretical < t_lo - margin || theoretical > t_hi + margin {
                    assert!(
                        observed < o_lo || observed > o_hi,
                        "{tol:?} t={theoretical}"
                    );
                }
            }
        }
    }
}
