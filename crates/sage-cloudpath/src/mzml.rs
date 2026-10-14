use crate::numpress::{NumpressCodec, NumpressError};
use flate2::read::ZlibDecoder;
use quick_xml::events::{BytesEnd, BytesStart, BytesText, Event};
use quick_xml::Reader;
use sage_core::spectrum::{Precursor, Representation};
use sage_core::{mass::Tolerance, spectrum::RawSpectrum};
use std::io::Read;
use tokio::io::AsyncBufRead;

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
/// Which tag are we inside?
enum State {
    Spectrum,
    Scan,
    BinaryDataArray,
    Binary,
    Precursor,
    SelectedIon,
}

#[derive(Copy, Clone, Debug)]
enum BinaryKind {
    Intensity,
    Mz,
    Noise,
}

#[derive(Copy, Clone, Debug)]
enum Dtype {
    F32,
    F64,
}

// MUST supply only one of the following
const ZLIB_COMPRESSION: &[u8] = b"MS:1000574";
const NO_COMPRESSION: &[u8] = b"MS:1000576";
const NUMPRESS_LINEAR: &[u8] = b"MS:1002312";
const NUMPRESS_PIC: &[u8] = b"MS:1002313";
const NUMPRESS_SLOF: &[u8] = b"MS:1002314";
const NUMPRESS_LINEAR_ZLIB: &[u8] = b"MS:1002746";
const NUMPRESS_PIC_ZLIB: &[u8] = b"MS:1002747";
const NUMPRESS_SLOF_ZLIB: &[u8] = b"MS:1002748";

// MUST supply only one of the following
const INTENSITY_ARRAY: &[u8] = b"MS:1000515";
const MZ_ARRAY: &[u8] = b"MS:1000514";
const NOISE_ARRAY: &[u8] = b"MS:1002744";

// MUST supply only one of the following
const FLOAT_64: &[u8] = b"MS:1000523";
const FLOAT_32: &[u8] = b"MS:1000521";

const MS_LEVEL: &[u8] = b"MS:1000511";
const PROFILE: &[u8] = b"MS:1000128";
const CENTROID: &[u8] = b"MS:1000127";
const TOTAL_ION_CURRENT: &[u8] = b"MS:1000285";

const SCAN_START_TIME: &[u8] = b"MS:1000016";
const UNIT_SECONDS: &[u8] = b"UO:0000010";
const UNIT_MINUTES: &[u8] = b"UO:0000031";
const ION_INJECTION_TIME: &[u8] = b"MS:1000927";

const SELECTED_ION_MZ: &[u8] = b"MS:1000744";
const SELECTED_ION_INT: &[u8] = b"MS:1000042";
const SELECTED_ION_CHARGE: &[u8] = b"MS:1000041";

const ISO_WINDOW_TARGET: &[u8] = b"MS:1000827";
const ISO_WINDOW_LOWER: &[u8] = b"MS:1000828";
const ISO_WINDOW_UPPER: &[u8] = b"MS:1000829";

const INVERSE_ION_MOBILITY: &[u8] = b"MS:1002815";

pub struct MzMLReader {
    ms_level: Option<u8>,
    // If set to Some(level) and noise intensities are present in the MzML file,
    // divide intensities at this MS-level by noise to calculate S/N
    signal_to_noise: Option<u8>,

    file_id: usize,
}

impl MzMLReader {
    /// Create a new [`MzMlReader`] with a minimum MS level filter
    ///
    /// # Example
    ///
    /// A minimum level of 2 will not parse or return MS1 scans
    pub fn with_file_id_and_level_filter(file_id: usize, ms_level: u8) -> Self {
        Self {
            ms_level: Some(ms_level),
            file_id,
            signal_to_noise: None,
        }
    }

    pub fn with_file_id(file_id: usize) -> Self {
        Self {
            ms_level: None,
            signal_to_noise: None,
            file_id,
        }
    }

    pub fn set_file_id(&mut self, file_id: usize) -> &mut Self {
        self.file_id = file_id;
        self
    }

    pub fn set_signal_to_noise(&mut self, sn: Option<u8>) -> &mut Self {
        self.signal_to_noise = sn;
        self
    }

    /// Every spectrum of a whole mzML document.
    pub async fn parse<B: AsyncBufRead + Unpin>(
        &self,
        b: B,
    ) -> Result<Vec<RawSpectrum>, MzMLError> {
        let mut reader = Reader::from_reader(b);
        let mut buf = Vec::new();
        let mut builder = SpectrumBuilder::new(self);
        let mut spectra = Vec::new();
        loop {
            match reader.read_event_into_async(&mut buf).await {
                Ok(Event::Eof) => break,
                Ok(event) => {
                    if let Some(spectrum) = builder.event(&event)? {
                        spectra.push(spectrum);
                    }
                }
                Err(err) => {
                    log::error!("unhandled XML error while parsing mzML: {}", err)
                }
            }
            buf.clear();
        }
        Ok(spectra)
    }

    /// The spectrum whose `<spectrum>` element starts `xml`, or `None` when this reader
    /// skips it (MS level filter, no ion current). Bytes after its `</spectrum>` are
    /// ignored, so `xml` may run up to the next spectrum's start.
    pub fn parse_spectrum(&self, xml: &[u8]) -> Result<Option<RawSpectrum>, MzMLError> {
        let mut reader = Reader::from_reader(xml);
        let mut buf = Vec::new();
        let mut builder = SpectrumBuilder::new(self);
        loop {
            match reader.read_event_into(&mut buf)? {
                Event::Eof => return Ok(None),
                Event::End(ref ev) if ev.name().into_inner() == b"spectrum" => {
                    return builder.end(ev);
                }
                event => {
                    builder.event(&event)?;
                }
            }
            buf.clear();
        }
    }
}

macro_rules! extract {
    ($ev:expr, $key:expr) => {
        $ev.try_get_attribute($key)?
            .ok_or(MzMLError::Malformed)?
            .value
    };
}

macro_rules! extract_value {
    ($ev:expr) => {{
        let s = $ev
            .try_get_attribute(b"value")?
            .ok_or(MzMLError::Malformed)?
            .value;
        std::str::from_utf8(&s)?.parse()?
    }};
}

/// Builds spectra from the XML events inside their `<spectrum>` elements. All state
/// resets at each `<spectrum>`, so a spectrum parses the same alone (`parse_spectrum`)
/// as within its whole file (`parse`).
///
/// Here be dragons - this is a giant imperative, stateful event handler.
struct SpectrumBuilder<'r> {
    reader: &'r MzMLReader,
    state: Option<State>,
    compression: bool,
    binary_dtype: Dtype,
    binary_array: Option<BinaryKind>,
    numpress: Option<NumpressCodec>,
    spectrum: RawSpectrum,
    precursor: Precursor,
    iso_window_lo: Option<f32>,
    iso_window_hi: Option<f32>,
    noise_array: Vec<f32>,
    inflated: Vec<u8>,
}

impl<'r> SpectrumBuilder<'r> {
    fn new(reader: &'r MzMLReader) -> Self {
        Self {
            reader,
            state: None,
            compression: false,
            binary_dtype: Dtype::F64,
            binary_array: None,
            numpress: None,
            spectrum: RawSpectrum::default_with_file_id(reader.file_id),
            precursor: Precursor::default(),
            iso_window_lo: None,
            iso_window_hi: None,
            noise_array: Vec::new(),
            inflated: Vec::with_capacity(4096),
        }
    }

    /// A finished spectrum when `event` closes one this reader keeps.
    fn event(&mut self, event: &Event) -> Result<Option<RawSpectrum>, MzMLError> {
        match event {
            Event::Start(ev) => self.start(ev)?,
            Event::Empty(ev) => self.empty(ev)?,
            Event::Text(text) => self.text(text)?,
            Event::End(ev) => return self.end(ev),
            _ => {}
        }
        Ok(None)
    }

    fn start(&mut self, ev: &BytesStart) -> Result<(), MzMLError> {
        // State transition into child tag
        self.state = match (ev.name().into_inner(), self.state) {
            (b"spectrum", _) => {
                *self = Self::new(self.reader);
                Some(State::Spectrum)
            }
            (b"scan", Some(State::Spectrum)) => Some(State::Scan),
            (b"binaryDataArray", Some(State::Spectrum)) => {
                self.numpress = None;
                Some(State::BinaryDataArray)
            }
            (b"binary", Some(State::BinaryDataArray)) => Some(State::Binary),
            (b"precursor", Some(State::Spectrum)) => Some(State::Precursor),
            (b"selectedIon", Some(State::Precursor)) => Some(State::SelectedIon),
            _ => self.state,
        };
        match ev.name().into_inner() {
            b"spectrum" => {
                let id = extract!(ev, b"id");
                let id = std::str::from_utf8(&id)?;
                self.spectrum.id = id.to_string();
            }
            b"precursor" => {
                // Not all precursor fields have a spectrumRef
                if let Some(scan) = ev.try_get_attribute(b"spectrumRef")? {
                    let scan = std::str::from_utf8(&scan.value)?;
                    self.precursor.spectrum_ref = Some(scan.to_string())
                }
            }
            _ => {}
        }
        Ok(())
    }

    fn empty(&mut self, ev: &BytesStart) -> Result<(), MzMLError> {
        match (self.state, ev.name().into_inner()) {
            (Some(State::BinaryDataArray), b"cvParam") => {
                let accession = extract!(ev, b"accession");
                match accession.as_ref() {
                    ZLIB_COMPRESSION => self.compression = true,
                    NO_COMPRESSION => self.compression = false,
                    NUMPRESS_LINEAR => self.numpress = Some(NumpressCodec::Linear),
                    NUMPRESS_PIC => self.numpress = Some(NumpressCodec::PositiveInteger),
                    NUMPRESS_SLOF => self.numpress = Some(NumpressCodec::ShortLoggedFloat),
                    NUMPRESS_LINEAR_ZLIB => {
                        self.numpress = Some(NumpressCodec::Linear);
                        self.compression = true;
                    }
                    NUMPRESS_PIC_ZLIB => {
                        self.numpress = Some(NumpressCodec::PositiveInteger);
                        self.compression = true;
                    }
                    NUMPRESS_SLOF_ZLIB => {
                        self.numpress = Some(NumpressCodec::ShortLoggedFloat);
                        self.compression = true;
                    }
                    FLOAT_64 => self.binary_dtype = Dtype::F64,
                    FLOAT_32 => self.binary_dtype = Dtype::F32,
                    INTENSITY_ARRAY => self.binary_array = Some(BinaryKind::Intensity),
                    MZ_ARRAY => self.binary_array = Some(BinaryKind::Mz),
                    NOISE_ARRAY => self.binary_array = Some(BinaryKind::Noise),
                    _ => {
                        // Unknown CV - perhaps noise
                        self.binary_array = None;
                    }
                }
            }
            (Some(State::Spectrum), b"cvParam") => {
                let accession = extract!(ev, b"accession");
                match accession.as_ref() {
                    MS_LEVEL => {
                        let level = extract_value!(ev);
                        if let Some(filter) = self.reader.ms_level {
                            if level != filter {
                                self.spectrum =
                                    RawSpectrum::default_with_file_id(self.reader.file_id);
                                self.state = None;
                            }
                        }
                        self.spectrum.ms_level = level;
                    }
                    PROFILE => self.spectrum.representation = Representation::Profile,
                    CENTROID => self.spectrum.representation = Representation::Centroid,
                    TOTAL_ION_CURRENT => {
                        let value = extract_value!(ev);
                        if value == 0.0 {
                            // No ion current, break out of current self.state
                            self.spectrum = RawSpectrum::default_with_file_id(self.reader.file_id);
                            self.state = None;
                        } else {
                            self.spectrum.total_ion_current = value;
                        }
                    }
                    _ => {}
                }
            }
            (Some(State::Precursor), b"cvParam") => {
                let accession = extract!(ev, b"accession");
                match accession.as_ref() {
                    ISO_WINDOW_TARGET => {
                        // use isolation window target for self.precursor m/z, e.g. to handle
                        // DIA setups where the mzML conversion software doesn't write
                        // a selection ion tag
                        if self.precursor.mz == 0.0 {
                            self.precursor.mz = extract_value!(ev)
                        }
                    }
                    ISO_WINDOW_LOWER => self.iso_window_lo = Some(extract_value!(ev)),
                    ISO_WINDOW_UPPER => self.iso_window_hi = Some(extract_value!(ev)),
                    _ => {}
                }
            }
            (Some(State::SelectedIon), b"cvParam") => {
                let accession = extract!(ev, b"accession");
                match accession.as_ref() {
                    SELECTED_ION_CHARGE => {
                        self.precursor.charge = Some(extract_value!(ev));
                    }
                    SELECTED_ION_MZ => {
                        let val = extract_value!(ev);
                        if val != 0.0 {
                            self.precursor.mz = val;
                        }
                    }
                    SELECTED_ION_INT => {
                        self.precursor.intensity = Some(extract_value!(ev));
                    }
                    INVERSE_ION_MOBILITY => {
                        self.precursor.inverse_ion_mobility = Some(extract_value!(ev));
                    }
                    _ => {}
                }
            }
            (Some(State::Scan), b"cvParam") => {
                let accession = extract!(ev, b"accession");
                match accession.as_ref() {
                    SCAN_START_TIME => {
                        let scan_start_time = extract_value!(ev);
                        let unit = extract!(ev, b"unitAccession");

                        self.spectrum.scan_start_time = match unit.as_ref() {
                            UNIT_SECONDS => scan_start_time / 60.0,
                            UNIT_MINUTES => scan_start_time,
                            _ => return Err(MzMLError::Malformed),
                        };
                    }
                    ION_INJECTION_TIME => {
                        self.spectrum.ion_injection_time = extract_value!(ev);
                    }
                    INVERSE_ION_MOBILITY => {
                        self.precursor.inverse_ion_mobility = Some(extract_value!(ev));
                    }
                    _ => {}
                }
            }

            _ => {}
        }
        Ok(())
    }

    fn text(&mut self, text: &BytesText) -> Result<(), MzMLError> {
        if let Some(State::Binary) = self.state {
            if let Some(filter) = self.reader.ms_level {
                if self.spectrum.ms_level != filter {
                    return Ok(());
                }
            }
            let raw = text.unescape()?;
            // There are occasionally empty binary data arrays, or unknown CVs
            if raw.is_empty() || self.binary_array.is_none() {
                return Ok(());
            }
            let decoded = base64::decode(raw.as_bytes())?;
            let bytes = match self.compression {
                false => &decoded,
                true => {
                    self.inflated.clear();
                    ZlibDecoder::new(decoded.as_slice()).read_to_end(&mut self.inflated)?;
                    &self.inflated
                }
            };

            let array = match (self.numpress, self.binary_dtype) {
                (Some(codec), _) => {
                    let mut values = Vec::new();
                    codec.decode(bytes, &mut values)?;
                    values
                }
                (None, Dtype::F32) => {
                    let mut buf: [u8; 4] = [0; 4];
                    bytes
                        .chunks(4)
                        .filter(|chunk| chunk.len() == 4)
                        .map(|chunk| {
                            buf.copy_from_slice(chunk);
                            f32::from_le_bytes(buf)
                        })
                        .collect::<Vec<f32>>()
                }
                (None, Dtype::F64) => {
                    let mut buf: [u8; 8] = [0; 8];
                    bytes
                        .chunks(8)
                        .map(|chunk| {
                            buf.copy_from_slice(chunk);
                            f64::from_le_bytes(buf) as f32
                        })
                        .collect::<Vec<f32>>()
                }
            };

            match self.binary_array {
                Some(BinaryKind::Intensity) => {
                    self.spectrum.intensity = array;
                }
                Some(BinaryKind::Mz) => {
                    self.spectrum.mz = array;
                }
                Some(BinaryKind::Noise) => {
                    self.noise_array = array;
                }
                None => {}
            }

            self.binary_array = None;
        }
        Ok(())
    }

    fn end(&mut self, ev: &BytesEnd) -> Result<Option<RawSpectrum>, MzMLError> {
        let mut finished = None;
        self.state = match (self.state, ev.name().into_inner()) {
            (Some(State::Binary), b"binary") => Some(State::BinaryDataArray),
            (Some(State::BinaryDataArray), b"binaryDataArray") => Some(State::Spectrum),
            (Some(State::SelectedIon), b"selectedIon") => Some(State::Precursor),
            (Some(State::Precursor), b"precursor") => {
                if self.precursor.mz != 0.0 {
                    let mut precursor = std::mem::take(&mut self.precursor);
                    precursor.isolation_window = match (self.iso_window_lo, self.iso_window_hi) {
                        // Offsets bound the theoretical m/z to
                        // [target − lower, target + upper], i.e.
                        // observed − theoretical ∈ [−upper, lower].
                        (Some(lower), Some(upper)) => Some(Tolerance::Da(-upper, lower)),
                        _ => None,
                    };
                    self.spectrum.precursors.push(precursor);
                }
                Some(State::Spectrum)
            }
            (Some(State::Scan), b"scan") => Some(State::Spectrum),
            (_, b"spectrum") => {
                let spectrum = std::mem::replace(
                    &mut self.spectrum,
                    RawSpectrum::default_with_file_id(self.reader.file_id),
                );
                let allow = self
                    .reader
                    .ms_level
                    .map(|level| level == spectrum.ms_level)
                    .unwrap_or(true);
                finished = match (allow, self.reader.signal_to_noise) {
                    (true, Some(level))
                        if level == spectrum.ms_level && !self.noise_array.is_empty() =>
                    {
                        let mut spectrum = spectrum;
                        spectrum
                            .intensity
                            .iter_mut()
                            .zip(self.noise_array.iter())
                            .for_each(|(int, noise)| *int /= noise);
                        self.noise_array.clear();
                        Some(spectrum)
                    }
                    (true, _) => Some(spectrum),
                    (false, _) => None,
                };
                None
            }
            _ => self.state,
        };
        Ok(finished)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum MzMLError {
    #[error("malformed MzML")]
    Malformed,
    #[error("unsupported cvParam {0}")]
    UnsupportedCV(String),
    #[error("XML parsing error: {0}")]
    XMLError(#[from] quick_xml::Error),
    #[error("io error: {0}")]
    IOError(#[from] std::io::Error),
    #[error("utf8 error: {0}")]
    Utf8Error(#[from] std::str::Utf8Error),
    #[error("error parsing float: {0}")]
    FloatError(#[from] std::num::ParseFloatError),
    #[error("error parsing int: {0}")]
    IntError(#[from] std::num::ParseIntError),
    #[error("error decoding base64: {0}")]
    Base64Error(#[from] base64::DecodeError),
    #[error("error decoding numpress: {0}")]
    NumpressError(#[from] NumpressError),
}

#[cfg(test)]
mod test {
    use sage_core::{mass::Tolerance, spectrum::Representation};

    use super::{MzMLError, MzMLReader};

    /// Real-data check of numpress decoding: `SAGE_NUMPRESS_PAIR=<numpress.mzML>,<zlib.mzML>`
    /// (the same spectra written both ways, e.g. by pmsms2mzml). Intensities must match
    /// exactly (pic stores integers), m/z within numpress' linear-prediction rounding.
    #[test]
    #[ignore]
    fn numpress_mzml_matches_its_zlib_twin() {
        let pair =
            std::env::var("SAGE_NUMPRESS_PAIR").expect("SAGE_NUMPRESS_PAIR=<numpress>,<zlib>");
        let (numpress, zlib) = pair.split_once(',').expect("two comma-separated paths");
        let numpress = crate::util::read_mzml(numpress, 0, None).expect("numpress mzML");
        let zlib = crate::util::read_mzml(zlib, 0, None).expect("zlib mzML");
        assert_eq!(numpress.len(), zlib.len());
        let mut max_ppm = 0f64;
        let mut peaks = 0usize;
        for (a, b) in numpress.iter().zip(&zlib) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.intensity, b.intensity, "intensities of {}", a.id);
            assert_eq!(a.mz.len(), b.mz.len());
            for (&x, &y) in a.mz.iter().zip(&b.mz) {
                max_ppm = max_ppm.max(((x as f64 - y as f64) / y as f64).abs() * 1e6);
            }
            peaks += a.mz.len();
        }
        eprintln!(
            "{} spectra, {peaks} peaks, max m/z difference {max_ppm:.6} ppm",
            numpress.len()
        );
        assert!(max_ppm < 0.01, "max m/z difference {max_ppm} ppm");
    }

    #[tokio::test]
    async fn parse_spectrum_issue_78() -> Result<(), MzMLError> {
        let s = r#"
        <spectrum id="spectrum=2442" index="286" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <cvParam cvRef="MS" accession="MS:1000528" name="lowest observed m/z" value="147.290603637695"/>
            <cvParam cvRef="MS" accession="MS:1000527" name="highest observed m/z" value="769.255798339843977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <userParam name="preset scan configuration" type="xsd:string" value="2"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                            <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();

        assert_eq!(s.id, "spectrum=2442");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert!(match s.precursors[0].inverse_ion_mobility {
            Some(x) => (x - 1.0786) < 0.0001,
            None => false,
        });
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-0.75, 1.5))
        );
        assert!((s.scan_start_time - 25.066).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        Ok(())
    }

    #[tokio::test]
    async fn parse_spectrum_issue_117() -> Result<(), MzMLError> {
        // The issue was that some converters write the ion mobility as part of the selected ion (as in the last test)
        // and some write it as part of the scan, as in this test. This test checks that it can be read
        // fom the scan section.
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <cvParam cvRef="MS" accession="MS:1000294" name="mass spectrum" />
            <cvParam cvRef="MS" accession="MS:1000130" name="positive scan" />
            <cvParam cvRef="MS" accession="MS:1000504" name="base peak m/z" value="638.352905273437955"/>
            <cvParam cvRef="MS" accession="MS:1000505" name="base peak intensity" value="113.885513305664006"/>
            <cvParam cvRef="MS" accession="MS:1000285" name="total ion current" value="793.395202636718977"/>
            <userParam name="filter string" type="xsd:string" value="ITMS + c NSI d w Full ms2 457.72@cid35.00 [115.00-930.00]"/>
            <scanList count="1">
                <cvParam cvRef="MS" accession="MS:1000795" name="no combination" />
                <scan >
                    <cvParam cvRef="MS" accession="MS:1000016" name="scan start time" value="1503.96166992188" unitAccession="UO:0000010" unitName="second" unitCvRef="UO" />
                    <cvParam cvRef="MS" accession="MS:1002815" name="inverse reduced ion mobility" value="1.078628" unitAccession="MS:1002814" unitName="volt-second per square centimeter"/>
                    <userParam name="[Thermo Trailer Extra]Monoisotopic M/Z:" type="xsd:double" value="457.723968505858977"/>
                    <scanWindowList count="1">
                        <scanWindow>
                            <cvParam cvRef="MS" accession="MS:1000501" name="scan window lower limit" value="115" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000500" name="scan window upper limit" value="930" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        </scanWindow>
                    </scanWindowList>
                </scan>
            </scanList>
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <activation>
                        <cvParam cvRef="MS" accession="MS:1000133" name="collision-induced dissociation" />
                        <cvParam cvRef="MS" accession="MS:1000045" name="collision energy" value="35.0"/>
                    </activation>
                </precursor>
            </precursorList>
            <binaryDataArrayList count="2">
                <binaryDataArray encodedLength="1088">
                    <cvParam cvRef="MS" accession="MS:1000514" name="m/z array" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    <cvParam cvRef="MS" accession="MS:1000523" name="64-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>AAAAoExpYkAAAACA3MpkQAAAAACph2VAAAAAAE4wZkAAAACAlMdmQAAAAECZAmdAAAAAwP9jaEAAAADgj4ZoQAAAAGC7HWlAAAAAAOXFaUAAAADg+4dqQAAAAMC1pmpAAAAA4IGFa0AAAACAaUZsQAAAACBzYW1AAAAAANCjbUAAAACAQ6duQAAAAIDsxG5AAAAAQKIlb0AAAACA5z9vQAAAAIDuw29AAAAAAJQicEAAAAAg9UZwQAAAAKCeVHBAAAAAIInEcEAAAACAcs5wQAAAAOA6BHFAAAAAADoOcUAAAAAgfcRxQAAAAOA68nFAAAAAoPExckAAAADATKVyQAAAAMC10nJAAAAAwBJHc0AAAAAA7FNzQAAAAIAYkXNAAAAAgJzRc0AAAABgE2R0QAAAAMCrc3RAAAAAgE+zdEAAAAAAhMR0QAAAAIC64XRAAAAA4Cf/dEAAAADgy3B1QAAAAMCVgnVAAAAAoDugdUAAAACAX/Z1QAAAAAAAB3ZAAAAAgO4XdkAAAABAqEJ2QAAAAIDp8nZAAAAAIAgRd0AAAACggzR3QAAAAODwT3dAAAAAIHJsd0AAAAAA4YJ3QAAAAGC91ndAAAAAAL3id0AAAADg0xZ4QAAAAOA5NXhAAAAAYDaPeEAAAACgK7p4QAAAACCm0XhAAAAA4GHkeEAAAADgyPJ4QAAAAOB5/3hAAAAAoFtNeUAAAADA8H15QAAAAGAHtXlAAAAAoD7HeUAAAAAAEtR5QAAAAGCx5XlAAAAA4NEJekAAAAAgtVN6QAAAACDCX3pAAAAAIAqmekAAAACg4OR6QAAAAGDymnxAAAAAICV/fUAAAAAgd6Z9QAAAAKDYA4BAAAAAoCoVgEAAAACA/kOAQAAAAKCpYoBAAAAA4MycgEAAAADA3DyBQAAAAKCbrIFAAAAAoPC6gUAAAADgV22CQAAAACABY4NAAAAAQE+qg0AAAADA0vKDQAAAAEDz+oNAAAAAoIxrhEAAAADg6euEQAAAAIAuDIVAAAAAoOwjhUAAAACgZUuFQAAAAADdm4VAAAAAoCzrh0AAAABgYvWHQAAAAOALCohA</binary>
                </binaryDataArray>
                <binaryDataArray encodedLength="544">
                    <cvParam cvRef="MS" accession="MS:1000515" name="intensity array" unitAccession="MS:1000131" unitName="number of detector counts" unitCvRef="MS"/>
                    <cvParam cvRef="MS" accession="MS:1000521" name="32-bit float" />
                    <cvParam cvRef="MS" accession="MS:1000576" name="no compression" />
                    <binary>3FlbQDg/ZUB8w3FAV2fMQMiOnkCXfP4/T2I2QC6qskAnhOZA/NU2QCc2QEAI1UhAQcAbQRrziUBmHq5AXutSQWZDbkAZGWdAzt6lQYNptUDSFDNBoY4IQAYaQEDeT7Q/16HGP9GtXUCITrQ/Rxu0Pzhc6j9mpjZAX1X8P7tPQ0AqxS5BZTzZPye+m0B7Sa5AfPsPQRr/W0CYwBRBwDh3QMAmtD/nq6E/bJHGPxJ9UUDsy/dAoCYMQRM2a0BkAR9Boo5pQMV0VEArYu5A4kaMQAyTI0BQPRJAML3TQCKVCED85+tArObGP1BVP0EtJuVAdyKAQFjctkFQa2NBixMTQXyyjUFX8eo/IHelQTdFcEFo1zZAhagsQAO53EBIugRB0M+gQfhBgkH0MsJAbGlIQZXg+EHe6CZBsbA2QHMHOECtW6BAjE2oQUpZckBasZ1AtKl3QEZYIUHkip1AQX7TQPqF60GNuaE/USk2QGLF40Im65ZAmXqlQBGuSUC70KBAAneMQeK3aEB87MVA5NigQE/Wb0BO475A</binary>
                </binaryDataArray>
            </binaryDataArrayList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!(match s.precursors[0].inverse_ion_mobility {
            Some(x) => (x - 1.0786) < 0.0001,
            None => false,
        });

        // The rest of these assertions just make sure the integrity of the spectrum is maintained
        assert_eq!(s.id, "spectrum=8678309");
        assert_eq!(s.ms_level, 2);
        assert_eq!(s.representation, Representation::Centroid);
        assert_eq!(s.precursors.len(), 1);
        assert_eq!(s.precursors[0].charge, Some(2));
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-0.75, 1.5))
        );
        assert!((s.scan_start_time - 25.066).abs() < 0.0001);
        assert_eq!(s.ion_injection_time, 0.0);
        assert_eq!(s.intensity.len(), s.mz.len());
        Ok(())
    }

    #[tokio::test]
    async fn spectra_parse_alone_as_in_their_file() -> Result<(), MzMLError> {
        let with_window = r#"<spectrum id="first" index="0" defaultArrayLength="0">
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1"><precursor>
                <isolationWindow>
                    <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" />
                    <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" />
                </isolationWindow>
                <selectedIonList count="1"><selectedIon>
                    <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.72" />
                </selectedIon></selectedIonList>
            </precursor></precursorList>
        </spectrum>
        "#;
        let without_window = r#"<spectrum id="second" index="1" defaultArrayLength="0">
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1"><precursor>
                <selectedIonList count="1"><selectedIon>
                    <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="512.3" />
                </selectedIon></selectedIonList>
            </precursor></precursorList>
        </spectrum>
        "#;
        let file = format!(
            "<mzML><run><spectrumList>{with_window}{without_window}</spectrumList></run></mzML>"
        );
        let reader = MzMLReader::with_file_id(0);
        let whole = reader.parse(file.as_bytes()).await?;

        assert_eq!(whole.len(), 2);
        assert_eq!(
            whole[0].precursors[0].isolation_window,
            Some(Tolerance::Da(-0.75, 1.5))
        );
        // The first spectrum's isolation window doesn't leak into the second.
        assert_eq!(whole[1].precursors[0].isolation_window, None);

        let second_start = file.find("<spectrum id=\"second\"").unwrap();
        let first_start = file.find("<spectrum id=\"first\"").unwrap();
        let first = reader.parse_spectrum(&file.as_bytes()[first_start..second_start])?;
        let second = reader.parse_spectrum(&file.as_bytes()[second_start..])?;
        assert_eq!(format!("{:?}", first.unwrap()), format!("{:?}", whole[0]));
        assert_eq!(format!("{:?}", second.unwrap()), format!("{:?}", whole[1]));
        Ok(())
    }

    #[tokio::test]
    async fn parse_spectrum_issue_210() -> Result<(), MzMLError> {
        // Handle cases where both isolation window target m/z is set and different than selected ion m/z
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-0.75, 1.5))
        );

        // Check different ordering of fields in mzML
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000744" name="selected ion m/z" value="457.723968505859" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.723968) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-0.75, 1.5))
        );

        // Check fallback keeping iso window m/z of fields in mzML
        let s = r#"
        <spectrum id="spectrum=8678309" index="8678309" defaultArrayLength="102" dataProcessingRef="dp_sp_1">
            <cvParam cvRef="MS" accession="MS:1000127" name="centroid spectrum" />
            <cvParam cvRef="MS" accession="MS:1000511" name="ms level" value="2" />
            <precursorList count="1">
                <precursor>
                    <isolationWindow>
                        <cvParam cvRef="MS" accession="MS:1000827" name="isolation window target m/z" value="457.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000828" name="isolation window lower offset" value="1.5" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                        <cvParam cvRef="MS" accession="MS:1000829" name="isolation window upper offset" value="0.75" unitAccession="MS:1000040" unitName="m/z" unitCvRef="MS" />
                    </isolationWindow>
                    <selectedIonList count="1">
                        <selectedIon>
                            <cvParam cvRef="MS" accession="MS:1000041" name="charge state" value="2" />
                        </selectedIon>
                    </selectedIonList>
                </precursor>
            </precursorList>
        </spectrum>
        "#;
        let mut spectra = MzMLReader::with_file_id(0).parse(s.as_bytes()).await?;

        assert_eq!(spectra.len(), 1);
        let s = spectra.pop().unwrap();
        assert!((s.precursors[0].mz - 457.75) < 0.0001);
        assert_eq!(
            s.precursors[0].isolation_window,
            Some(Tolerance::Da(-0.75, 1.5))
        );
        Ok(())
    }
}
