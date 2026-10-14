//! Random access into a gzip file: zlib's zran method (`examples/zran.c`). One full
//! decompression records access points every `spacing` decompressed bytes; any offset is
//! then reached from the last point before it by inflating at most `spacing` bytes.
//! An access point is a clone of miniz_oxide's inflater (its bit buffer included) plus
//! the last 32 KB of output, which stands in for zlib's `inflatePrime` and
//! `inflateSetDictionary`. Members of multi-member files (`pigz`, `bgzip`) are followed
//! transparently, and the indexing pass checks every member's CRC-32 and length.

use miniz_oxide::inflate::core::{decompress, DecompressorOxide};
use miniz_oxide::inflate::TINFLStatus;

/// DEFLATE's largest back-reference distance; the output ring holds that much history.
const RING: usize = 1 << 15;
/// Decompressed bytes between access points: ~43 KB per point, at most this many bytes
/// re-inflated to reach an offset.
pub const DEFAULT_SPACING: u64 = 4 << 20;

#[derive(Debug, thiserror::Error)]
pub enum GzipError {
    #[error("not a gzip member at compressed byte {0}")]
    NotGzip(usize),
    #[error("gzip data ends inside a member")]
    Truncated,
    #[error("corrupt deflate data: {0:?}")]
    Corrupt(TINFLStatus),
    #[error("gzip member ending at compressed byte {0} fails its CRC-32 or length check")]
    Checksum(usize),
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Phase {
    Header,
    Inflating,
    Trailer,
}

/// A position in the decompressed stream, resumable on its own.
#[derive(Clone)]
struct Inflater {
    phase: Phase,
    state: Box<DecompressorOxide>,
    ring: Box<[u8]>,
    ring_pos: usize,
    /// Compressed bytes consumed.
    input: usize,
    /// Decompressed bytes produced.
    output: u64,
}

enum Step<'r> {
    Data(&'r [u8]),
    MemberEnd { crc: u32, size: u32, end: usize },
    End,
}

impl Inflater {
    fn new() -> Self {
        Self {
            phase: Phase::Header,
            state: Box::default(),
            ring: vec![0; RING].into_boxed_slice(),
            ring_pos: 0,
            input: 0,
            output: 0,
        }
    }

    /// The next bytes of output (up to the ring's end), or the end of a member or file.
    fn step<'r>(&'r mut self, data: &[u8]) -> Result<Step<'r>, GzipError> {
        match self.phase {
            Phase::Header => {
                if self.input >= data.len() {
                    return Ok(Step::End);
                }
                self.input += member_header_len(data, self.input)?;
                *self.state = DecompressorOxide::new();
                self.phase = Phase::Inflating;
                Ok(Step::Data(&[]))
            }
            Phase::Trailer => {
                let trailer = data
                    .get(self.input..self.input + 8)
                    .ok_or(GzipError::Truncated)?;
                let crc = u32::from_le_bytes(trailer[..4].try_into().expect("4 bytes"));
                let size = u32::from_le_bytes(trailer[4..].try_into().expect("4 bytes"));
                self.input += 8;
                self.phase = Phase::Header;
                Ok(Step::MemberEnd {
                    crc,
                    size,
                    end: self.input,
                })
            }
            Phase::Inflating => {
                let start = self.ring_pos;
                let (status, consumed, written) = decompress(
                    &mut self.state,
                    &data[self.input..],
                    &mut self.ring,
                    start,
                    0,
                );
                self.input += consumed;
                self.output += written as u64;
                self.ring_pos = (start + written) & (RING - 1);
                match status {
                    TINFLStatus::Done => self.phase = Phase::Trailer,
                    TINFLStatus::HasMoreOutput => {}
                    TINFLStatus::NeedsMoreInput => return Err(GzipError::Truncated),
                    other => return Err(GzipError::Corrupt(other)),
                }
                Ok(Step::Data(&self.ring[start..start + written]))
            }
        }
    }
}

impl Inflater {
    /// The first output offset still held by the ring (it holds the last `RING` bytes).
    fn retained_from(&self) -> u64 {
        self.output.saturating_sub(RING as u64)
    }

    /// Appends output `start..end`, which the ring still holds, to `out`.
    fn copy_retained(&self, start: u64, end: u64, out: &mut Vec<u8>) {
        let back = (self.output - start) as usize;
        let from = (self.ring_pos + RING - back) % RING;
        let len = (end - start) as usize;
        let first = len.min(RING - from);
        out.extend_from_slice(&self.ring[from..from + first]);
        out.extend_from_slice(&self.ring[..len - first]);
    }
}

/// Length of the gzip member header (RFC 1952) at `at`.
fn member_header_len(data: &[u8], at: usize) -> Result<usize, GzipError> {
    const FHCRC: u8 = 2;
    const FEXTRA: u8 = 4;
    const FNAME: u8 = 8;
    const FCOMMENT: u8 = 16;
    let header = &data[at..];
    if header.len() < 10 || header[..3] != [0x1f, 0x8b, 8] {
        return Err(GzipError::NotGzip(at));
    }
    let flags = header[3];
    let mut len = 10;
    if flags & FEXTRA != 0 {
        let extra = header.get(len..len + 2).ok_or(GzipError::Truncated)?;
        len += 2 + u16::from_le_bytes([extra[0], extra[1]]) as usize;
    }
    for field in [FNAME, FCOMMENT] {
        if flags & field != 0 {
            let terminator = header
                .get(len..)
                .and_then(|rest| rest.iter().position(|&b| b == 0))
                .ok_or(GzipError::Truncated)?;
            len += terminator + 1;
        }
    }
    if flags & FHCRC != 0 {
        len += 2;
    }
    if len > header.len() {
        return Err(GzipError::Truncated);
    }
    Ok(len)
}

/// Access points into one gzip file's decompressed stream.
pub struct GzipIndex {
    points: Vec<Inflater>,
    len: u64,
    spacing: u64,
}

impl GzipIndex {
    /// Decompresses `data` once, handing every decompressed chunk with its offset to
    /// `on_data`, and records an access point every `spacing` decompressed bytes.
    pub fn build(
        data: &[u8],
        spacing: u64,
        mut on_data: impl FnMut(u64, &[u8]),
    ) -> Result<Self, GzipError> {
        let mut inflater = Inflater::new();
        let mut points = vec![inflater.clone()];
        let mut crc = crc32fast::Hasher::new();
        let mut member_size = 0u32;
        loop {
            let offset = inflater.output;
            match inflater.step(data)? {
                Step::Data(bytes) => {
                    crc.update(bytes);
                    member_size = member_size.wrapping_add(bytes.len() as u32);
                    on_data(offset, bytes);
                }
                Step::MemberEnd {
                    crc: stored,
                    size,
                    end,
                } => {
                    let computed = std::mem::replace(&mut crc, crc32fast::Hasher::new());
                    if computed.finalize() != stored || member_size != size {
                        return Err(GzipError::Checksum(end));
                    }
                    member_size = 0;
                }
                Step::End => break,
            }
            let last = points.last().expect("starts with one point").output;
            if inflater.output - last >= spacing {
                points.push(inflater.clone());
            }
        }
        Ok(Self {
            points,
            len: inflater.output,
            spacing,
        })
    }

    /// Decompressed length.
    pub fn len(&self) -> u64 {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// A reader of `data` (the file this index was built from), positioned nowhere yet.
    pub fn cursor<'a>(&'a self, data: &'a [u8]) -> GzipCursor<'a> {
        GzipCursor {
            index: self,
            data,
            inflater: None,
            jumps: 0,
        }
    }

    fn point_before(&self, offset: u64) -> &Inflater {
        let after = self.points.partition_point(|point| point.output <= offset);
        &self.points[after - 1]
    }
}

/// Reads decompressed byte ranges, continuing from its last read when the range starts
/// in the bytes it still holds or less than one spacing ahead.
pub struct GzipCursor<'a> {
    index: &'a GzipIndex,
    data: &'a [u8],
    inflater: Option<Inflater>,
    /// Reads that restarted from an access point.
    jumps: usize,
}

impl GzipCursor<'_> {
    /// Replaces `out` with decompressed bytes `start..end`.
    pub fn read(&mut self, start: u64, end: u64, out: &mut Vec<u8>) -> Result<(), GzipError> {
        out.clear();
        let index = self.index;
        let inflater = match self.inflater.take() {
            Some(inflater)
                if inflater.retained_from() <= start && start < inflater.output + index.spacing =>
            {
                inflater
            }
            _ => {
                self.jumps += 1;
                index.point_before(start).clone()
            }
        };
        let inflater = self.inflater.insert(inflater);
        if start < inflater.output {
            inflater.copy_retained(start, end.min(inflater.output), out);
        }
        while inflater.output < end {
            let offset = inflater.output;
            match inflater.step(self.data)? {
                Step::Data(bytes) => {
                    let chunk_end = offset + bytes.len() as u64;
                    if chunk_end > start {
                        let from = start.saturating_sub(offset) as usize;
                        let to = (end.min(chunk_end) - offset) as usize;
                        out.extend_from_slice(&bytes[from..to]);
                    }
                }
                Step::MemberEnd { .. } => {}
                Step::End => break,
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn text(lines: usize) -> Vec<u8> {
        (0..lines)
            .map(|i| format!("line {i}: {}\n", "abcdefgh".repeat(i % 13)))
            .collect::<String>()
            .into_bytes()
    }

    fn gzip(members: &[&[u8]]) -> Vec<u8> {
        let mut out = Vec::new();
        for member in members {
            let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::default());
            encoder.write_all(member).unwrap();
            out.extend(encoder.finish().unwrap());
        }
        out
    }

    #[test]
    fn reads_any_range_of_single_and_multi_member_files() {
        let plain = text(60_000);
        let (first, second) = plain.split_at(plain.len() / 3);
        for compressed in [gzip(&[&plain]), gzip(&[first, second])] {
            let mut streamed = Vec::new();
            let index = GzipIndex::build(&compressed, 64 * 1024, |offset, bytes| {
                assert_eq!(offset, streamed.len() as u64);
                streamed.extend_from_slice(bytes);
            })
            .unwrap();
            assert_eq!(streamed, plain);
            assert_eq!(index.len(), plain.len() as u64);
            assert!(index.points.len() > 10);

            let mut cursor = index.cursor(&compressed);
            let mut out = Vec::new();
            let n = plain.len() as u64;
            let mut state = 12345u64;
            let mut ranges: Vec<(u64, u64)> = (0..200)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    let start = (state >> 11) % n;
                    (start, (start + 40_000).min(n))
                })
                .collect();
            ranges.extend([(0, 10), (n - 10, n), (n / 3 - 5, n / 3 + 5)]);
            for (start, end) in ranges {
                cursor.read(start, end, &mut out).unwrap();
                assert_eq!(out, &plain[start as usize..end as usize], "{start}..{end}");
            }
        }
    }

    #[test]
    fn consecutive_reads_continue_without_jumping() {
        let plain = text(60_000);
        let (first, second) = plain.split_at(plain.len() / 3);
        let compressed = gzip(&[first, second]);
        let index = GzipIndex::build(&compressed, 64 * 1024, |_, _| {}).unwrap();
        let mut cursor = index.cursor(&compressed);
        let mut out = Vec::new();
        let n = plain.len() as u64;
        let mut start = 0;
        while start < n {
            let end = (start + 7_919).min(n);
            cursor.read(start, end, &mut out).unwrap();
            assert_eq!(out, &plain[start as usize..end as usize]);
            start = end;
        }
        assert_eq!(cursor.jumps, 1);
    }

    #[test]
    fn a_corrupted_trailer_fails_the_check() {
        let mut compressed = gzip(&[&text(1000)]);
        let crc_at = compressed.len() - 8;
        compressed[crc_at] ^= 1;
        assert!(matches!(
            GzipIndex::build(&compressed, DEFAULT_SPACING, |_, _| {}),
            Err(GzipError::Checksum(_))
        ));
    }

    #[test]
    fn a_truncated_file_is_an_error() {
        let compressed = gzip(&[&text(1000)]);
        assert!(GzipIndex::build(
            &compressed[..compressed.len() / 2],
            DEFAULT_SPACING,
            |_, _| {}
        )
        .is_err());
    }
}
