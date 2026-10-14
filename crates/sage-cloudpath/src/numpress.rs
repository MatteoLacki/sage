//! MS-Numpress decoding (linear prediction, positive integer, short logged float),
//! ported from the reference implementation (ms-numpress `MSNumpress.cpp`) with every
//! read bounds-checked. Values are computed in f64 as there, then stored as f32.

#[derive(Debug, PartialEq, Eq, thiserror::Error)]
pub enum NumpressError {
    #[error("numpress array of {0} bytes is too short")]
    TooShort(usize),
    #[error("corrupt numpress array")]
    Corrupt,
}

/// Which numpress codec an mzML binary array uses.
#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum NumpressCodec {
    Linear,
    PositiveInteger,
    ShortLoggedFloat,
}

impl NumpressCodec {
    /// Appends the values decoded from `data` to `out`.
    pub fn decode(self, data: &[u8], out: &mut Vec<f32>) -> Result<(), NumpressError> {
        match self {
            Self::Linear => decode_linear(data, out),
            Self::PositiveInteger => decode_pic(data, out),
            Self::ShortLoggedFloat => decode_slof(data, out),
        }
    }
}

const FIXED_POINT_BYTES: usize = 8;

/// The fixed point is stored big-endian.
fn fixed_point(data: &[u8]) -> Result<f64, NumpressError> {
    let bytes = data
        .get(..FIXED_POINT_BYTES)
        .ok_or(NumpressError::TooShort(data.len()))?;
    Ok(f64::from_be_bytes(bytes.try_into().expect("8 bytes")))
}

/// Integers stored as variable runs of half-bytes, high half of each byte first.
struct HalfBytes<'a> {
    data: &'a [u8],
    byte: usize,
    in_low_half: bool,
}

impl<'a> HalfBytes<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            byte: 0,
            in_low_half: false,
        }
    }

    /// No integer left: past the last byte, or at a last low half-byte of zero padding.
    fn at_end(&self) -> bool {
        let last = self.data.len().saturating_sub(1);
        self.byte >= self.data.len()
            || (self.byte == last && self.in_low_half && self.data[last] & 0xf == 0)
    }

    fn next_half_byte(&mut self) -> Result<u8, NumpressError> {
        let byte = *self.data.get(self.byte).ok_or(NumpressError::Corrupt)?;
        if self.in_low_half {
            self.byte += 1;
            self.in_low_half = false;
            Ok(byte & 0xf)
        } else {
            self.in_low_half = true;
            Ok(byte >> 4)
        }
    }

    /// A head half-byte `h` <= 8 says `h` leading zero half-bytes were dropped, `h` > 8
    /// that `h - 8` leading 0xf half-bytes were; the remaining half-bytes follow, lowest
    /// first.
    fn next_int(&mut self) -> Result<u32, NumpressError> {
        let head = self.next_half_byte()? as usize;
        let (dropped, mut value) = match head {
            0..=8 => (head, 0u32),
            _ => {
                let ones = head - 8;
                let leading_ones = (0..ones).fold(0u32, |v, i| v | (0xf000_0000 >> (4 * i)));
                (ones, leading_ones)
            }
        };
        if dropped == 8 {
            return Ok(value);
        }
        let pending_low_half = usize::from(!self.in_low_half);
        if self.byte + ((8 - dropped) - pending_low_half) / 2 >= self.data.len() {
            return Err(NumpressError::Corrupt);
        }
        for i in dropped..8 {
            let half_byte = self.next_half_byte()? as u32;
            value |= half_byte << ((i - dropped) * 4);
        }
        Ok(value)
    }
}

/// Linear prediction (m/z, retention time): a fixed point, the first two values as
/// little-endian u32, then each value's difference from `2·prev − prev_prev` as an int.
pub fn decode_linear(data: &[u8], out: &mut Vec<f32>) -> Result<(), NumpressError> {
    if data.len() == FIXED_POINT_BYTES {
        return Ok(());
    }
    let fixed_point = fixed_point(data)?;
    let seed = |at: usize| -> Result<i64, NumpressError> {
        let bytes = data
            .get(at..at + 4)
            .ok_or(NumpressError::TooShort(data.len()))?;
        Ok(u32::from_le_bytes(bytes.try_into().expect("4 bytes")) as i64)
    };
    let mut before_previous: i64;
    let mut previous = seed(8)?;
    out.push((previous as f64 / fixed_point) as f32);
    if data.len() == 12 {
        return Ok(());
    }
    let mut current = seed(12)?;
    out.push((current as f64 / fixed_point) as f32);

    let mut ints = HalfBytes::new(&data[16..]);
    while !ints.at_end() {
        before_previous = previous;
        previous = current;
        let difference = ints.next_int()? as i32;
        current = previous
            .checked_sub(before_previous)
            .and_then(|step| previous.checked_add(step))
            .and_then(|extrapolated| extrapolated.checked_add(difference as i64))
            .ok_or(NumpressError::Corrupt)?;
        out.push((current as f64 / fixed_point) as f32);
    }
    Ok(())
}

/// Positive integer (intensities): one half-byte int per value.
pub fn decode_pic(data: &[u8], out: &mut Vec<f32>) -> Result<(), NumpressError> {
    let mut ints = HalfBytes::new(data);
    while !ints.at_end() {
        out.push(ints.next_int()? as f32);
    }
    Ok(())
}

/// Short logged float (intensities): a fixed point, then little-endian u16 values `x`
/// standing for `exp(x / fixed_point) − 1`.
pub fn decode_slof(data: &[u8], out: &mut Vec<f32>) -> Result<(), NumpressError> {
    let fixed_point = fixed_point(data)?;
    let values = &data[FIXED_POINT_BYTES..];
    if values.len() % 2 != 0 {
        return Err(NumpressError::Corrupt);
    }
    out.extend(values.chunks_exact(2).map(|pair| {
        let x = u16::from_le_bytes([pair[0], pair[1]]);
        ((x as f64 / fixed_point).exp() - 1.0) as f32
    }));
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;

    // Encoded and decoded by the reference C++ (ms-numpress MSNumpress.cpp, vendored in
    // git/pmsms2mzml); decoded values as f32 bits.
    const LINEAR: &[(&str, &[u8], &[u32])] = &[
        (
            "linear_0",
            &[0x40, 0xc3, 0x88, 0x00, 0x00, 0x00, 0x00, 0x00],
            &[],
        ),
        (
            "linear_1",
            &[
                0x41, 0x52, 0x67, 0x64, 0xc0, 0x00, 0x00, 0x00, 0x03, 0xff, 0xff, 0x7f,
            ],
            &[0x43de8fcd],
        ),
        (
            "linear_2",
            &[
                0x41, 0x52, 0x58, 0xde, 0x40, 0x00, 0x00, 0x00, 0xa5, 0xf8, 0x9a, 0x7f, 0x8b, 0xfe,
                0xff, 0x7f,
            ],
            &[0x43de8fcd, 0x43df4000],
        ),
        (
            "linear_3",
            &[
                0x41, 0x39, 0x2f, 0xd0, 0x00, 0x00, 0x00, 0x00, 0x40, 0xad, 0xd6, 0x09, 0xf8, 0x74,
                0xfc, 0x09, 0x34, 0xfb, 0x46, 0x1a, 0x24, 0xeb, 0x59, 0x05, 0x1f, 0xab, 0x56, 0x70,
                0xab, 0x95, 0xc1, 0x08, 0x00, 0xc8, 0x36, 0xc3, 0x10,
            ],
            &[
                0x42c80000, 0x42cb0000, 0x42ce8000, 0x43482000, 0x44bb8001, 0x44bb8002, 0x44d49f9b,
            ],
        ),
        (
            "linear_4",
            &[
                0x41, 0x42, 0x34, 0x56, 0x00, 0x00, 0x00, 0x00, 0xa6, 0xc7, 0x75, 0x1c, 0xdb, 0x24,
                0x76, 0x1c, 0x88, 0x1c, 0x61, 0xa5, 0x4e, 0x95, 0xb5, 0x3c, 0xc1, 0x09, 0xd9, 0x72,
                0xdf, 0x70, 0x3d, 0xee, 0x15, 0x08, 0x09, 0x2a, 0xa9, 0x8f, 0x30,
            ],
            &[
                0x43481c29, 0x43481eb8, 0x43482148, 0x434823d7, 0x43964000, 0x43968000, 0x44960800,
                0x44963800, 0x44ce4000,
            ],
        ),
    ];
    const PIC: &[(&str, &[u8], &[u32])] = &[
        ("pic_0", &[], &[]),
        ("pic_1", &[0x80], &[0x00000000]),
        (
            "pic_2",
            &[
                0x71, 0x7f, 0x60, 0x16, 0xff, 0x50, 0x01, 0x4f, 0xff, 0xf3, 0x00, 0x00, 0x10, 0x00,
                0x49, 0x53, 0x77,
            ],
            &[
                0x3f800000, 0x41700000, 0x41800000, 0x437f0000, 0x43800000, 0x477fff00, 0x47800000,
                0x4eee6b28,
            ],
        ),
        (
            "pic_3",
            &[
                0x7c, 0x87, 0x36, 0xd4, 0x50, 0x04, 0x50, 0x08, 0x75, 0x76, 0x77,
            ],
            &[
                0x41400000, 0x00000000, 0x40400000, 0x429a0000, 0x44800000, 0x45000000, 0x40a00000,
                0x40c00000, 0x40e00000,
            ],
        ),
    ];
    const SLOF: &[(&str, &[u8], &[u32])] = &[
        (
            "slof_0",
            &[0x40, 0x8f, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00],
            &[],
        ),
        (
            "slof_1",
            &[
                0x40, 0xb6, 0x3c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x60, 0x14, 0x9d, 0x66,
                0x79, 0xd1, 0xfc, 0xff,
            ],
            &[0x00000000, 0x3fc006cf, 0x42c7fd95, 0x4640e57b, 0x47c35343],
        ),
        (
            "slof_2",
            &[
                0x40, 0xb2, 0x87, 0x00, 0x00, 0x00, 0x00, 0x00, 0xaf, 0x19, 0xaf, 0x19, 0xaf, 0x19,
                0xf7, 0xff,
            ],
            &[0x403ffd51, 0x403ffd51, 0x403ffd51, 0x49742461],
        ),
    ];

    fn cases() -> impl Iterator<Item = (NumpressCodec, &'static str, &'static [u8], &'static [u32])>
    {
        let tag = |codec, table: &'static [(&'static str, &'static [u8], &'static [u32])]| {
            table
                .iter()
                .map(move |&(name, bytes, bits)| (codec, name, bytes, bits))
        };
        tag(NumpressCodec::Linear, LINEAR)
            .chain(tag(NumpressCodec::PositiveInteger, PIC))
            .chain(tag(NumpressCodec::ShortLoggedFloat, SLOF))
    }

    fn decoded_bits(codec: NumpressCodec, bytes: &[u8]) -> Result<Vec<u32>, NumpressError> {
        let mut out = Vec::new();
        codec.decode(bytes, &mut out)?;
        Ok(out.iter().map(|v| v.to_bits()).collect())
    }

    #[test]
    fn decodes_exactly_as_the_reference() {
        for (codec, name, bytes, expected) in cases() {
            assert_eq!(
                decoded_bits(codec, bytes).as_deref(),
                Ok(expected),
                "{name}"
            );
        }
    }

    #[test]
    fn truncated_arrays_fail_or_decode_a_prefix() {
        for (codec, name, bytes, expected) in cases() {
            for length in 0..bytes.len() {
                if let Ok(bits) = decoded_bits(codec, &bytes[..length]) {
                    assert!(expected.starts_with(&bits), "{name} cut to {length} bytes");
                }
            }
        }
    }

    #[test]
    fn garbage_never_panics() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next_byte = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 24) as u8
        };
        for length in 0..200 {
            let bytes: Vec<u8> = (0..length).map(|_| next_byte()).collect();
            for codec in [
                NumpressCodec::Linear,
                NumpressCodec::PositiveInteger,
                NumpressCodec::ShortLoggedFloat,
            ] {
                let _ = decoded_bits(codec, &bytes);
            }
        }
    }

    #[test]
    fn odd_short_logged_float_payload_is_corrupt() {
        let (_, _, bytes, _) = cases().find(|case| case.1 == "slof_1").unwrap();
        assert_eq!(
            decoded_bits(NumpressCodec::ShortLoggedFloat, &bytes[..bytes.len() - 1]),
            Err(NumpressError::Corrupt)
        );
    }
}
