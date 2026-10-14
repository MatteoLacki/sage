//! mzML and MGF files read one spectrum at a time through an index of where each
//! spectrum's text starts: found by a parallel byte search of a plain file, or during
//! the one decompression pass of a gzipped one (`gzip_index`), which also records the
//! access points later reads start from.

use crate::gzip_index::{GzipCursor, GzipError, GzipIndex, DEFAULT_SPACING};
use crate::mgf::{DefaultParams, MgfReader};
use crate::mzml::{MzMLError, MzMLReader};
use memmap2::Mmap;
use rayon::prelude::*;
use sage_core::spectrum::RawSpectrum;
use std::fs::File;
use std::path::Path;

/// MGF parameters set before the first spectrum are read from at most this much text.
const MAX_HEADER_BYTES: usize = 1 << 20;
/// Bytes kept before the scan position of a decompressing scan, for MGF's
/// start-of-line check.
const LOOKBEHIND: usize = 256;

#[derive(Debug, thiserror::Error)]
pub enum TextSpectraError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("gzip error: {0}")]
    Gzip(#[from] GzipError),
    #[error("mzML error: {0}")]
    MzML(#[from] MzMLError),
    #[error("MGF spectrum is not UTF-8: {0}")]
    Utf8(#[from] std::str::Utf8Error),
}

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum TextFormat {
    MzML,
    Mgf,
}

impl TextFormat {
    /// The text every spectrum starts with.
    fn opening(self) -> &'static [u8] {
        match self {
            Self::MzML => b"<spectrum",
            Self::Mgf => b"BEGIN IONS",
        }
    }

    /// Whether the opening at `at` in `text` starts a spectrum, or `None` if `text` ends
    /// too early to tell. `text[0]` is the file's first byte when `at_file_start`.
    fn starts_spectrum(self, text: &[u8], at: usize, at_file_start: bool) -> Option<bool> {
        match self {
            // `<spectrum id=…`, not `<spectrumList`.
            Self::MzML => {
                let next = *text.get(at + self.opening().len())?;
                Some(next.is_ascii_whitespace() || next == b'>')
            }
            // `BEGIN IONS` first on its line, after optional blanks.
            Self::Mgf => {
                let before = &text[..at];
                match before.iter().rposition(|&b| b == b'\n') {
                    Some(newline) => Some(
                        before[newline + 1..]
                            .iter()
                            .all(|&b| b == b' ' || b == b'\t'),
                    ),
                    None => Some(at_file_start && before.iter().all(|&b| b == b' ' || b == b'\t')),
                }
            }
        }
    }
}

/// Finds spectrum starts in decompressed text arriving in chunks.
struct StreamScan {
    format: TextFormat,
    pending: Vec<u8>,
    /// Decompressed offset of `pending[0]`.
    pending_offset: u64,
    /// Openings starting before this offset are decided.
    decided: u64,
    starts: Vec<u64>,
    header: Vec<u8>,
}

impl StreamScan {
    fn new(format: TextFormat) -> Self {
        Self {
            format,
            pending: Vec::new(),
            pending_offset: 0,
            decided: 0,
            starts: Vec::new(),
            header: Vec::new(),
        }
    }

    fn feed(&mut self, bytes: &[u8]) {
        if self.starts.is_empty() && self.header.len() < MAX_HEADER_BYTES {
            self.header.extend_from_slice(bytes);
        }
        self.pending.extend_from_slice(bytes);
        self.scan(false);
    }

    fn finish(mut self) -> (Vec<u64>, Vec<u8>) {
        self.scan(true);
        if let Some(&first) = self.starts.first() {
            self.header.truncate(first as usize);
        }
        (self.starts, self.header)
    }

    fn scan(&mut self, at_end: bool) {
        let opening = self.format.opening();
        let from = (self.decided - self.pending_offset) as usize;
        let mut undecided = None;
        for found in memchr::memmem::find_iter(&self.pending[from..], opening) {
            let at = from + found;
            match self
                .format
                .starts_spectrum(&self.pending, at, self.pending_offset == 0)
            {
                Some(true) => self.starts.push(self.pending_offset + at as u64),
                Some(false) => {}
                None if at_end => {}
                None => {
                    undecided = Some(at);
                    break;
                }
            }
        }
        let searched_to = self.pending.len().saturating_sub(opening.len() - 1);
        let decided_to = undecided.unwrap_or(searched_to).max(from);
        self.decided = self.pending_offset + decided_to as u64;
        let keep_from = decided_to.saturating_sub(LOOKBEHIND);
        self.pending.drain(..keep_from);
        self.pending_offset += keep_from as u64;
    }
}

/// Spectrum starts in a whole plain file, searched in parallel ranges.
fn scan_plain(text: &[u8], format: TextFormat) -> Vec<u64> {
    const RANGE: usize = 8 << 20;
    let opening = format.opening();
    let ranges = text.len().div_ceil(RANGE);
    (0..ranges)
        .into_par_iter()
        .flat_map_iter(|range| {
            let from = range * RANGE;
            let to = (from + RANGE).min(text.len());
            let search_to = (to + opening.len() - 1).min(text.len());
            memchr::memmem::find_iter(&text[from..search_to], opening)
                .map(move |found| from + found)
                .filter(move |&at| at < to)
                .filter(move |&at| format.starts_spectrum(text, at, true).unwrap_or(false))
                .map(|at| at as u64)
        })
        .collect()
}

/// An mzML or MGF file (plain or gzipped) whose spectra are read by index.
pub struct TextSpectra {
    format: TextFormat,
    file: Mmap,
    gzip: Option<GzipIndex>,
    starts: Vec<u64>,
    /// Offset just past the last spectrum's text (the end of the decompressed file).
    end: u64,
    mzml: MzMLReader,
    mgf: MgfReader,
    mgf_defaults: DefaultParams,
}

/// One thread's reading state: a reused text buffer and, for gzip, a decompressor that
/// moves forward from one spectrum to the next.
pub struct TextReader<'a> {
    text: Vec<u8>,
    cursor: Option<GzipCursor<'a>>,
}

impl TextSpectra {
    /// Indexes `path`'s MS2 spectra; MS1 and higher levels are skipped when read.
    pub fn open(path: &Path, file_id: usize, format: TextFormat) -> Result<Self, TextSpectraError> {
        // SAFETY: the file is only read, and must not change while this run reads it.
        let file = unsafe { Mmap::map(&File::open(path)?)? };
        let is_gzip = file.starts_with(&[0x1f, 0x8b]);
        let (gzip, starts, header, end) = if is_gzip {
            let mut scan = StreamScan::new(format);
            let index = GzipIndex::build(&file, DEFAULT_SPACING, |_, bytes| scan.feed(bytes))?;
            let (starts, header) = scan.finish();
            let end = index.len();
            (Some(index), starts, header, end)
        } else {
            let starts = scan_plain(&file, format);
            let header_end = starts.first().map_or(file.len(), |&s| s as usize);
            let header = file[..header_end.min(MAX_HEADER_BYTES)].to_vec();
            (None, starts, header, file.len() as u64)
        };
        let mgf = MgfReader::with_file_id(file_id);
        let mgf_defaults = match format {
            TextFormat::Mgf => mgf.parse_header(&String::from_utf8_lossy(&header)),
            TextFormat::MzML => DefaultParams::default_with_file_id(file_id),
        };
        Ok(Self {
            format,
            file,
            gzip,
            starts,
            end,
            mzml: MzMLReader::with_file_id_and_level_filter(file_id, 2),
            mgf,
            mgf_defaults,
        })
    }

    pub fn len(&self) -> usize {
        self.starts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    pub fn reader(&self) -> TextReader<'_> {
        TextReader {
            text: Vec::new(),
            cursor: None,
        }
    }

    /// Spectrum `index`, or `None` when its parser skips it.
    pub fn read<'a>(
        &'a self,
        index: usize,
        reader: &mut TextReader<'a>,
    ) -> Result<Option<RawSpectrum>, TextSpectraError> {
        let start = self.starts[index];
        let end = self.starts.get(index + 1).copied().unwrap_or(self.end);
        let text: &[u8] = match &self.gzip {
            None => &self.file[start as usize..end as usize],
            Some(index) => {
                let cursor = reader
                    .cursor
                    .get_or_insert_with(|| index.cursor(&self.file));
                cursor.read(start, end, &mut reader.text)?;
                &reader.text
            }
        };
        match self.format {
            TextFormat::MzML => Ok(self.mzml.parse_spectrum(text)?),
            TextFormat::Mgf => Ok(self
                .mgf
                .parse_block(&self.mgf_defaults, std::str::from_utf8(text)?)),
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use flate2::write::GzEncoder;
    use std::io::Write;

    fn mgf(spectra: usize) -> String {
        let mut text = String::from("# comment\nCHARGE=3+\n");
        for i in 0..spectra {
            text.push_str(&format!(
                "BEGIN IONS\nTITLE=s{i}\nPEPMASS={}\n{} 10\n{} 20\nEND IONS\n\n",
                400.0 + i as f32,
                100.0 + i as f32,
                200.0 + i as f32
            ));
        }
        text
    }

    fn write_file(dir: &Path, name: &str, bytes: &[u8]) -> std::path::PathBuf {
        let path = dir.join(name);
        std::fs::write(&path, bytes).unwrap();
        path
    }

    fn gzipped(bytes: &[u8]) -> Vec<u8> {
        let mut encoder = GzEncoder::new(Vec::new(), flate2::Compression::default());
        encoder.write_all(bytes).unwrap();
        encoder.finish().unwrap()
    }

    #[test]
    fn mgf_read_by_index_equals_whole_file_parse_plain_and_gzipped() {
        let dir = std::env::temp_dir().join(format!("sage_text_spectra_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let text = mgf(3000);
        let whole = MgfReader::with_file_id(7).parse(text.clone()).unwrap();
        for path in [
            write_file(&dir, "a.mgf", text.as_bytes()),
            write_file(&dir, "a.mgf.gz", &gzipped(text.as_bytes())),
        ] {
            let spectra = TextSpectra::open(&path, 7, TextFormat::Mgf).unwrap();
            assert_eq!(spectra.len(), whole.len());
            let mut reader = spectra.reader();
            for (index, expected) in whole.iter().enumerate() {
                let read = spectra.read(index, &mut reader).unwrap().unwrap();
                assert_eq!(format!("{read:?}"), format!("{expected:?}"));
            }
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn stream_scan_finds_openings_split_across_chunks() {
        let text = mgf(50);
        let expected = scan_plain(text.as_bytes(), TextFormat::Mgf);
        for chunk in [1, 3, 7, 64, 1000] {
            let mut scan = StreamScan::new(TextFormat::Mgf);
            text.as_bytes()
                .chunks(chunk)
                .for_each(|bytes| scan.feed(bytes));
            assert_eq!(scan.finish().0, expected, "chunks of {chunk}");
        }
        let mzml = "<mzML><spectrumList count=\"2\"><spectrum id=\"a\"></spectrum>\n<spectrum\tid=\"b\"></spectrum></spectrumList></mzML>";
        let expected = scan_plain(mzml.as_bytes(), TextFormat::MzML);
        assert_eq!(expected.len(), 2);
        for chunk in [1, 2, 5, 9, 10, 11] {
            let mut scan = StreamScan::new(TextFormat::MzML);
            mzml.as_bytes()
                .chunks(chunk)
                .for_each(|bytes| scan.feed(bytes));
            assert_eq!(scan.finish().0, expected, "chunks of {chunk}");
        }
    }
}
