//! One input file's MS2 spectra read one at a time by index, whatever their format, so
//! a search never holds a whole file's spectra (necromerge2 `plans/streaming_spectra.md`,
//! `docs/ai/streaming_spectra.md`).

use sage_cloudpath::pmsms::PmsmsSpectra;
use sage_cloudpath::tdf::TdfSpectra;
use sage_cloudpath::text_spectra::{TextReader, TextSpectra};
use sage_core::spectrum::{Precursor, ProcessedSpectrum, SpectrumProcessor};
use std::borrow::Cow;

pub enum IndexedSpectra {
    Pmsms(PmsmsSpectra),
    Tdf(TdfSpectra),
    /// mzML or MGF, plain or gzipped.
    Text(TextSpectra),
    /// Spectra a batch reader already read whole (TMT, LFQ and remote files).
    InMemory(Vec<ProcessedSpectrum>),
    /// A file that failed to open; logged, and searched as empty like the batch reader does.
    Unreadable,
}

/// One thread's state for reading `IndexedSpectra`: text files keep a buffer and a gzip
/// decompressor that moves forward through the thread's run of indices.
pub type IndexedReader<'a> = Option<TextReader<'a>>;

/// What planning a sliced search needs of a spectrum: its MS level and precursors.
pub struct SpectrumHeader {
    pub level: u8,
    pub precursors: Vec<Precursor>,
}

impl IndexedSpectra {
    pub fn len(&self) -> usize {
        match self {
            Self::Pmsms(spectra) => spectra.len(),
            Self::Tdf(spectra) => spectra.len(),
            Self::Text(spectra) => spectra.len(),
            Self::InMemory(spectra) => spectra.len(),
            Self::Unreadable => 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn reader(&self) -> IndexedReader<'_> {
        match self {
            Self::Text(spectra) => Some(spectra.reader()),
            _ => None,
        }
    }

    /// Spectrum `index`'s header, without processing its peaks (pmsms: without reading
    /// them); `None` for a spectrum the batch reader skips as well.
    pub fn header<'a>(
        &'a self,
        index: usize,
        reader: &mut IndexedReader<'a>,
    ) -> Option<SpectrumHeader> {
        let header = |level, precursors| SpectrumHeader { level, precursors };
        match self {
            Self::Pmsms(spectra) => Some(header(2, spectra.precursors_of(index))),
            Self::Tdf(spectra) => spectra
                .read(index)
                .map(|raw| header(raw.ms_level, raw.precursors)),
            Self::Text(spectra) => {
                let reader = reader
                    .as_mut()
                    .expect("text files read through a TextReader");
                match spectra.read(index, reader) {
                    Ok(raw) => raw.map(|raw| header(raw.ms_level, raw.precursors)),
                    Err(error) => {
                        log::error!("cannot read spectrum {index}: {error}");
                        None
                    }
                }
            }
            Self::InMemory(spectra) => Some(header(
                spectra[index].level,
                spectra[index].precursors.clone(),
            )),
            Self::Unreadable => None,
        }
    }

    /// Spectrum `index`, processed; `None` for one the batch reader skips as well.
    pub fn processed<'a>(
        &'a self,
        index: usize,
        processor: &SpectrumProcessor,
        reader: &mut IndexedReader<'a>,
    ) -> Option<Cow<'a, ProcessedSpectrum>> {
        match self {
            Self::Pmsms(spectra) => match spectra.read(index, |raw| processor.process(raw)) {
                Ok(spectrum) => Some(Cow::Owned(spectrum)),
                Err(error) => panic!("cannot read pmsms spectrum {index}: {error}"),
            },
            Self::Tdf(spectra) => spectra
                .read(index)
                .map(|raw| Cow::Owned(processor.process(raw))),
            Self::Text(spectra) => {
                let reader = reader
                    .as_mut()
                    .expect("text files read through a TextReader");
                match spectra.read(index, reader) {
                    Ok(raw) => raw.map(|raw| Cow::Owned(processor.process(raw))),
                    Err(error) => {
                        log::error!("cannot read spectrum {index}: {error}");
                        None
                    }
                }
            }
            Self::InMemory(spectra) => Some(Cow::Borrowed(&spectra[index])),
            Self::Unreadable => None,
        }
    }
}
