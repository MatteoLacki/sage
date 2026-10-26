//! Explicitly versioned, lossy base-peak-normalized fragment intensity codec.
//! Support/annotation IDs are separate: a zero code remains a retained channel.
use half::f16;

pub const SQRT8_BASEPEAK: &str = "sqrt8_basepeak_v1";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntensityEncoding {
    Fp16,
    Sqrt8Basepeak,
}

impl IntensityEncoding {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fp16 => "fp16",
            Self::Sqrt8Basepeak => SQRT8_BASEPEAK,
        }
    }
}

/// Encode a finite nonnegative value using its spectrum's maximum. The maximum
/// is transient and is deliberately not stored: decoded intensities are relative.
pub fn encode_sqrt8(value: f32, maximum: f32) -> Result<u8, &'static str> {
    if !value.is_finite() || !maximum.is_finite() || value < 0.0 || maximum < 0.0 || value > maximum
    {
        return Err("sqrt8 requires finite 0 <= intensity <= maximum");
    }
    if maximum == 0.0 {
        return Ok(0);
    }
    let scaled = (value / maximum).sqrt() * 255.0_f32;
    // Explicit nearest-even rounding, matching numpy.rint in FP32.
    let lower = scaled.floor();
    let remainder = scaled - lower;
    let mut code = lower as u8;
    if remainder > 0.5 || (remainder == 0.5 && code & 1 != 0) {
        code += 1;
    }
    Ok(code)
}

pub fn decode_sqrt8(code: u8) -> f32 {
    let root = code as f32 / 255.0_f32;
    root * root
}

/// Borrow encoded or legacy values without decoding the entire mapped library.
#[derive(Clone, Copy, Debug)]
pub enum FragmentIntensities<'a> {
    Fp16(&'a [f16]),
    Sqrt8Basepeak(&'a [u8]),
}

impl<'a> FragmentIntensities<'a> {
    pub fn len(self) -> usize {
        match self {
            Self::Fp16(values) => values.len(),
            Self::Sqrt8Basepeak(values) => values.len(),
        }
    }
    pub fn is_empty(self) -> bool {
        self.len() == 0
    }
    pub fn value(self, index: usize) -> f32 {
        match self {
            Self::Fp16(values) => values[index].to_f32(),
            Self::Sqrt8Basepeak(values) => decode_sqrt8(values[index]),
        }
    }
    pub fn fp16(self) -> Option<&'a [f16]> {
        match self {
            Self::Fp16(values) => Some(values),
            Self::Sqrt8Basepeak(_) => None,
        }
    }
    pub fn codes(self) -> Option<&'a [u8]> {
        match self {
            Self::Fp16(_) => None,
            Self::Sqrt8Basepeak(values) => Some(values),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn codebook_zero_range_and_scale_invariance() {
        for q in 0..=255 {
            assert_eq!(encode_sqrt8(decode_sqrt8(q), 1.0).unwrap(), q);
        }
        assert_eq!(encode_sqrt8(0.0, 0.0), Ok(0));
        assert_eq!(encode_sqrt8(-0.0, 0.0), Ok(0));
        assert_eq!(encode_sqrt8(8.9296875, 8.9296875), Ok(255));
        assert_eq!(encode_sqrt8(0.25, 1.0), Ok(128));
        assert_eq!(encode_sqrt8(2.0, 8.0), Ok(128));
        for (v, m) in [
            (f32::NAN, 1.0),
            (1.0, f32::INFINITY),
            (-1.0, 1.0),
            (2.0, 1.0),
            (1.0, 0.0),
        ] {
            assert!(encode_sqrt8(v, m).is_err());
        }
    }

    #[test]
    fn nearest_even_boundaries_and_normalized_error_bound() {
        for n in [0_u8, 1, 126, 127, 128, 253, 254] {
            let root = (n as f32 + 0.5) / 255.0;
            let value = root * root;
            let scaled = value.sqrt() * 255.0;
            if scaled == n as f32 + 0.5 {
                assert_eq!(
                    encode_sqrt8(value, 1.0).unwrap(),
                    if n & 1 == 0 { n } else { n + 1 }
                );
            }
        }
        for i in 0..=100000 {
            let value = i as f32 / 100000.0;
            let decoded = decode_sqrt8(encode_sqrt8(value, 1.0).unwrap());
            let bound = value.sqrt() / 255.0 + 1.0 / (4.0 * 255.0 * 255.0);
            assert!((decoded - value).abs() <= bound + 2e-7);
        }
    }
}
