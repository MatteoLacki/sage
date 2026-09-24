use crate::database::binary_search_slice;
use crate::mass::{Tolerance, NEUTRON, PROTON};

/// A charge-less peak at monoisotopic mass
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct Peak {
    pub intensity: f32,
    pub mass: f32,
}

/// Matched peak columns. Both arrays always have the same length and order.
#[derive(Clone, Default, Debug, PartialEq)]
pub struct PeakColumns {
    masses: Vec<f32>,
    intensities: Vec<f32>,
}

impl PeakColumns {
    pub fn len(&self) -> usize {
        self.masses.len()
    }
    pub fn masses(&self) -> &[f32] {
        &self.masses
    }
    pub fn intensities(&self) -> &[f32] {
        &self.intensities
    }
    pub fn peak(&self, index: usize) -> Peak {
        Peak {
            mass: self.masses[index],
            intensity: self.intensities[index],
        }
    }
    pub fn iter(&self) -> impl ExactSizeIterator<Item = Peak> + '_ {
        self.masses
            .iter()
            .zip(&self.intensities)
            .map(|(&mass, &intensity)| Peak { mass, intensity })
    }
}

impl From<Vec<Peak>> for PeakColumns {
    fn from(peaks: Vec<Peak>) -> Self {
        // The exact-size iterator lets unzip reserve the retained count in each column.
        let (masses, intensities) = peaks.into_iter().map(|p| (p.mass, p.intensity)).unzip();
        Self {
            masses,
            intensities,
        }
    }
}

impl Eq for Peak {}

impl PartialOrd for Peak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Peak {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.intensity
            .total_cmp(&other.intensity)
            .then_with(|| self.mass.total_cmp(&other.mass))
    }
}

/// A charge-less peak at monoisotopic mass with ion mobility
#[derive(PartialEq, Copy, Clone, Default, Debug)]
pub struct IMPeak {
    pub intensity: f32,
    pub mass: f32,
    pub mobility: f32, // I would use f16 but its not stable
}

impl Eq for IMPeak {}

impl PartialOrd for IMPeak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for IMPeak {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.intensity
            .total_cmp(&other.intensity)
            .then_with(|| self.mass.total_cmp(&other.mass))
            .then_with(|| self.mobility.total_cmp(&other.mobility))
    }
}

/// A de-isotoped peak, that might have some charge state information
#[derive(PartialEq, PartialOrd, Debug, Copy, Clone)]
pub struct Deisotoped {
    pub mz: f32,
    // Cumulative intensity of all isotopic peaks in the envelope higher than this one
    pub intensity: f32,
    // Assigned charge
    pub charge: Option<u8>,
    // If `Some(idx)`, idx is the index of the parent isotopic envelope
    pub envelope: Option<usize>,
}

#[derive(Debug, Clone)]
pub struct SpectrumProcessor {
    pub take_top_n: usize,
    pub min_deisotope_mz: f32,
    pub deisotope: bool,
    /// If `true`, trust that `RawSpectrum::mz` (and therefore the `Peak`s
    /// built from it) already arrives mass-ascending sorted, and replace
    /// `process()`'s own sort with a linear sortedness check
    /// (`is_sorted_by_mass`) that panics if that trust turns out misplaced.
    /// See `Search::assume_sorted_peaks` in `sage-cli` for the full
    /// rationale. Only affects paths where `process()` would otherwise sort
    /// at all -- see `is_sorted_by_mass`'s call site.
    pub assume_sorted_peaks: bool,
}

#[derive(Default, Debug, Clone)]
pub struct Precursor {
    pub mz: f32,
    pub intensity: Option<f32>,
    pub charge: Option<u8>,
    // pub scan: Option<usize>,
    pub spectrum_ref: Option<String>,
    pub isolation_window: Option<Tolerance>,
    pub inverse_ion_mobility: Option<f32>,
}

impl Precursor {
    /// This precursor's own tolerance if it has one, otherwise `default`.
    ///
    /// Lets individual precursors override the run-global precursor
    /// tolerance (e.g. per-precursor ppm bounds sourced from the pmsms
    /// precursors table) while falling back to existing behavior when
    /// `isolation_window` is unset.
    pub fn effective_precursor_tol(&self, default: Tolerance) -> Tolerance {
        self.isolation_window.unwrap_or(default)
    }
}

#[derive(Clone, Default, Debug)]
pub struct ProcessedSpectrum<T = PeakColumns> {
    /// MSn level
    pub level: u8,
    /// Scan ID
    pub id: String,
    /// File ID
    pub file_id: usize,
    /// Retention time in minutes
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Selected ions for precursors, if `level > 1`
    pub precursors: Vec<Precursor>,
    /// MS peaks, sorted by mass in ascending order
    pub peaks: T,
    /// Parallel to `peaks`. Non-empty only when deisotoping is enabled.
    /// `peak_charges[i]` is `peaks[i]`'s isotope-envelope-resolved charge,
    /// or `0` if deisotoping found no envelope for it (charge unknown, not
    /// "assumed 1" -- `peaks[i].mass` itself still uses scale 1 as its best
    /// available base in that case, see `SpectrumProcessor::process_ms2`;
    /// only this field's `0` records that the charge is actually unresolved,
    /// distinct from a genuinely confirmed charge-1 envelope, which is
    /// reported as `1`). Empty (all charges unknown) when deisotope=false.
    /// See `docs/ai/fragment_charge_hypothesis.md`.
    pub peak_charges: Vec<u8>,
    /// Total ion current
    pub total_ion_current: f32,
}

#[derive(Default, Debug, Clone)]
/// An unprocessed mass spectrum, as returned by a parser
/// *CRITICAL*: Users must set all fields manually, including `file_id`
pub struct RawSpectrum<M = Vec<f32>, I = Vec<f32>> {
    pub file_id: usize,
    /// MSn level
    pub ms_level: u8,
    /// Spectrum identifier
    pub id: String,
    /// Vector of precursors associated with this spectrum
    pub precursors: Vec<Precursor>,
    /// Profile or Centroided data
    pub representation: Representation,
    /// Scan start time in minutes
    pub scan_start_time: f32,
    /// Ion injection time
    pub ion_injection_time: f32,
    /// Total ion current
    pub total_ion_current: f32,
    /// M/z array
    pub mz: M,
    /// Intensity array
    pub intensity: I,
    /// Mobility array
    pub mobility: Option<Vec<f32>>,
}

impl RawSpectrum {
    /// Return a [`RawSpectrum`] with default values, but with the `file_id` field
    /// properly set
    pub fn default_with_file_id(file_id: usize) -> Self {
        Self {
            file_id,
            ..Default::default()
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
pub enum Representation {
    #[default]
    Profile,
    Centroid,
}

#[derive(Default)]
pub enum MS1Spectra {
    NoMobility(Vec<ProcessedSpectrum>),
    WithMobility(Vec<ProcessedSpectrum<Vec<IMPeak>>>),
    #[default]
    Empty,
}

/// The two peaks `select_matched_peaks` reports from one tolerance-window scan: the most
/// intense peak (what Sage actually matches a theoretical fragment to), and the peak closest
/// in mass to `center` (which may be a different, less intense peak). Named fields instead of
/// a `(usize, usize)` tuple so call sites read `m.most_intense`/`m.closest` rather than `.0`/`.1`.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PeakMatch {
    pub most_intense: usize,
    pub closest: usize,
}

/// Binary search followed by linear search to select, within `tolerance` of `center`, both the
/// most intense peak and the peak closest in mass to `center` -- a single scan, since both are
/// answered by the same tolerance window.
/// * `offset` - this parameter allows for a static adjustment to the lower and upper bounds of the search window.
///
/// Tie-breaks: most-intense keeps the highest-index peak on an exact intensity tie. Closest
/// prefers the smaller distance; on an exact distance tie, the more intense peak wins; on a tie
/// of both, the highest-index peak wins -- consistent with the most-intense tie-break.
///
/// Sage subtracts a proton (and assumes z=1) for all experimental peaks, and stores all fragments as monoisotopic
/// masses. This simplifies downstream calculations at multiple charge states, but it also subtly changes tolerance
/// bounds. For most applications this is completely OK to ignore - however, for exact similarity of TMT reporter ion
/// measurements with ProteomeDiscoverer, FragPipe, etc, we need to account for this minor difference (which has an impact
/// perhaps 0.01% of the time)
pub fn select_matched_peaks(
    peaks: &PeakColumns,
    center: f32,
    tolerance: Tolerance,
    offset: Option<f32>,
) -> Option<PeakMatch> {
    let (lo, hi) = tolerance.bounds(center);
    let (lo, hi) = (
        lo + offset.unwrap_or_default(),
        hi + offset.unwrap_or_default(),
    );

    let (i, j) = binary_search_slice(peaks.masses(), |mass, query| mass.total_cmp(query), lo, hi);

    let mut most_intense_idx = None;
    let mut max_int = 0.0;
    let mut closest_idx = None;
    let mut min_dist = f32::INFINITY;
    let mut closest_int = 0.0;
    // `binary_search_slice`'s left bound can overshoot by one element --
    // needed for correctness when it's used on a coarse per-bucket summary
    // (e.g. `page_search`'s outer search over `min_value`), where an exact
    // lower bound could skip a bucket that still straddles it. `peaks` is a
    // flat, fully-sorted array, so this re-check is redundant here, but
    // re-checking unconditionally keeps `binary_search_slice` one shared,
    // uniformly-safe helper instead of two near-identical variants.
    for idx in i..j {
        let mass = peaks.masses[idx];
        if mass < lo || mass > hi {
            continue;
        }
        let intensity = peaks.intensities[idx];

        if intensity >= max_int {
            max_int = intensity;
            most_intense_idx = Some(idx);
        }

        let dist = (mass - center).abs();
        if dist < min_dist || (dist == min_dist && intensity >= closest_int) {
            min_dist = dist;
            closest_int = intensity;
            closest_idx = Some(idx);
        }
    }

    Some(PeakMatch {
        most_intense: most_intense_idx?,
        closest: closest_idx?,
    })
}

// pub fn find_spectrum_by_id(
//     spectra: &[ProcessedSpectrum],
//     scan_id: usize,
// ) -> Option<&ProcessedSpectrum> {
//     // First try indexing by scan
//     if let Some(first) = spectra.get(scan_id.saturating_sub(1)) {
//         if first.scan == scan_id {
//             return Some(first);
//         }
//     }
//     // Fall back to binary search
//     let idx = spectra
//         .binary_search_by(|spec| spec.scan.cmp(&scan_id))
//         .ok()?;
//     spectra.get(idx)
// }

/// Raw intensities are converted to f32 before comparisons or accumulation.
pub trait RawIntensity: Copy {
    fn to_f32(self) -> f32;
}

impl RawIntensity for f32 {
    fn to_f32(self) -> f32 {
        self
    }
}

impl RawIntensity for u32 {
    fn to_f32(self) -> f32 {
        self as f32
    }
}

/// Deisotope a set of peaks by attempting to find C13 peaks under a given `ppm` tolerance
pub fn deisotope<T: RawIntensity>(
    mz: &[f32],
    int: &[T],
    max_charge: u8,
    ppm: f32,
    min_mz: f32,
) -> Vec<Deisotoped> {
    let mut peaks = mz
        .iter()
        .zip(int.iter())
        .map(|(mz, int)| Deisotoped {
            mz: *mz,
            intensity: int.to_f32(),
            envelope: None,
            charge: None,
        })
        .collect::<Vec<_>>();

    // Is the peak at index `i` an isotopic peak?
    for i in (0..mz.len()).rev() {
        // Two pointer approach, j is fast pointer
        let mut j = i.saturating_sub(1);
        while mz[i] - mz[j] <= NEUTRON + Tolerance::ppm_to_delta_mass(mz[i], ppm) && mz[j] >= min_mz
        {
            let delta = mz[i] - mz[j];
            let tol = Tolerance::ppm_to_delta_mass(mz[i], ppm);
            for charge in 1..=max_charge {
                let iso = NEUTRON / charge as f32;
                if (delta - iso).abs() <= tol && int[i].to_f32() < int[j].to_f32() {
                    // Make sure this peak isn't already part of an isotopic envelope
                    if let Some(existing) = peaks[i].charge {
                        if existing != charge {
                            continue;
                        }
                    }
                    peaks[j].intensity += peaks[i].intensity;
                    peaks[j].charge = Some(charge);
                    peaks[i].charge = Some(charge);
                    peaks[i].envelope = Some(j);
                }
            }
            j = j.saturating_sub(1);
            if j == 0 {
                break;
            }
        }
    }
    peaks
}

/// Path compression of isotopic envelope links
pub fn path_compression(peaks: &mut [Deisotoped]) {
    for idx in 0..peaks.len() {
        if let Some(parent) = peaks[idx].envelope {
            if let Some(upper) = peaks[parent].envelope {
                peaks[idx].envelope = Some(upper);
            }
            peaks[idx].intensity = 0.0;
        }
    }
}

impl ProcessedSpectrum {
    /// Compact masses, intensities and observed charges using the same mask.
    pub fn retain_peaks(&mut self, mut keep: impl FnMut(Peak) -> bool) {
        let mut dst = 0;
        for src in 0..self.peaks.len() {
            if keep(self.peaks.peak(src)) {
                self.peaks.masses[dst] = self.peaks.masses[src];
                self.peaks.intensities[dst] = self.peaks.intensities[src];
                if !self.peak_charges.is_empty() {
                    self.peak_charges[dst] = self.peak_charges[src];
                }
                dst += 1;
            }
        }
        self.peaks.masses.truncate(dst);
        self.peaks.intensities.truncate(dst);
        self.peak_charges.truncate(dst);
    }
}

impl<T> ProcessedSpectrum<T> {
    /// Was this spectrum's `peaks` built with isotope/charge deconvolution
    /// (`deisotope=true`)? Derived from `peak_charges` (empty iff
    /// `deisotope=false`, see that field's doc comment) rather than stored
    /// separately, so it can never drift out of sync with the peaks it
    /// describes.
    ///
    /// This distinction matters beyond just charge labeling: under
    /// `deisotope=true`, `process_ms2` merges each resolved isotope
    /// satellite's intensity into its monoisotopic root and drops the
    /// satellite from `peaks` entirely (`.filter(|peak|
    /// peak.envelope.is_none())`) -- so code that expects to find individual
    /// isotope-satellite peaks (e.g. `Scorer::observed_isotope_ladder`) only
    /// sees real data when this returns `false`.
    pub fn is_deisotoped(&self) -> bool {
        !self.peak_charges.is_empty()
    }

    pub fn extract_ms1_precursor(&self) -> Option<(f32, u8)> {
        let precursor = self.precursors.first()?;
        let charge = precursor.charge?;
        let mass = (precursor.mz - PROTON) * charge as f32;
        Some((mass, charge))
    }
    pub fn in_isolation_window(&self, mz: f32) -> Option<bool> {
        let precursor = self.precursors.first()?;
        let (lo, hi) = precursor.isolation_window?.bounds(precursor.mz - PROTON);
        Some(mz >= lo && mz <= hi)
    }
}

/// Is `peaks` already sorted ascending by `mass`? Hard precondition for
/// `select_matched_peaks`'s binary search, which returns silently wrong
/// or missing matches -- not a panic -- on unsorted input.
fn is_sorted_by_mass(peaks: &[Peak]) -> bool {
    peaks.windows(2).all(|w| w[0].mass <= w[1].mass)
}

impl SpectrumProcessor {
    /// Create a new [`SpectrumProcessor`]
    ///
    /// # Arguments
    /// * `take_top_n`: Keep only the top N most intense peaks from the spectrum
    /// * `min_fragment_mz`: Keep only fragments >= this m/z
    /// * `max_fragment_mz`: Keep only fragments <= this m/z
    /// * `deisotope`: Perform deisotoping & charge state deconvolution
    /// * `assume_sorted_peaks`: skip `process()`'s own mass-sort and instead
    ///   verify it with a linear scan, panicking if violated -- see the
    ///   field's own doc comment
    pub fn new(
        take_top_n: usize,
        deisotope: bool,
        min_deisotope_mz: f32,
        assume_sorted_peaks: bool,
    ) -> Self {
        Self {
            take_top_n,
            min_deisotope_mz,
            deisotope,
            assume_sorted_peaks,
        }
    }

    fn process_ms2<T: RawIntensity>(
        &self,
        should_deisotope: bool,
        spectrum: &RawSpectrum<impl AsRef<[f32]>, impl AsRef<[T]>>,
    ) -> (Vec<Peak>, Vec<u8>) {
        if spectrum.representation != Representation::Centroid {
            // Panic, because there's really nothing we can do with profile data
            panic!(
                "Scan {} contains profile data! Please convert to centroid",
                spectrum.id
            );
        }

        // If there is no precursor charge from the mzML file, then deisotope fragments up to z=3
        let charge = spectrum
            .precursors
            .first()
            .and_then(|p| p.charge)
            .unwrap_or(3);

        if should_deisotope {
            let mut peaks = deisotope(
                spectrum.mz.as_ref(),
                spectrum.intensity.as_ref(),
                charge,
                10.0,
                self.min_deisotope_mz,
            );
            peaks.sort_unstable_by(|a, b| {
                b.intensity
                    .total_cmp(&a.intensity)
                    .then_with(|| a.mz.total_cmp(&b.mz))
            });

            // Collect as (Peak, charge) pairs so we can sort both by mass together
            let mut pairs: Vec<(Peak, u8)> = peaks
                .into_iter()
                .filter(|peak| peak.envelope.is_none())
                .map(|peak| {
                    // Convert from MH* to M. `mass`'s scale always defaults
                    // to 1 when no envelope was resolved (the best available
                    // base for scoring.rs's fragment charge-hypothesis
                    // loop) -- but the charge *reported* below is 0 in that
                    // case, not 1, so it stays distinguishable from a
                    // genuinely confirmed charge-1 envelope. See
                    // docs/ai/fragment_charge_hypothesis.md.
                    let resolved_charge = peak.charge;
                    let scale_charge = resolved_charge.unwrap_or(1);
                    let mass = (peak.mz - PROTON) * scale_charge as f32;
                    (
                        Peak {
                            mass,
                            intensity: peak.intensity,
                        },
                        resolved_charge.unwrap_or(0),
                    )
                })
                .take(self.take_top_n)
                .collect();
            pairs.sort_by(|(a, _), (b, _)| a.mass.total_cmp(&b.mass));
            pairs.into_iter().unzip()
        } else {
            let mut peaks = spectrum
                .mz
                .as_ref()
                .iter()
                .zip(spectrum.intensity.as_ref().iter())
                .map(|(mz, &intensity)| {
                    let mass = (mz - PROTON) * 1.0;
                    Peak {
                        mass,
                        intensity: intensity.to_f32(),
                    }
                })
                .collect::<Vec<_>>();
            crate::heap::bounded_min_heapify(&mut peaks, self.take_top_n);
            peaks.truncate(self.take_top_n);
            (peaks, vec![])
        }
    }

    pub fn process<T: RawIntensity>(
        &self,
        spectrum: RawSpectrum<impl AsRef<[f32]>, impl AsRef<[T]>>,
    ) -> ProcessedSpectrum {
        let (mut peaks, peak_charges) = match spectrum.ms_level {
            2 => self.process_ms2(self.deisotope, &spectrum),
            _ => {
                let peaks = spectrum
                    .mz
                    .as_ref()
                    .iter()
                    .zip(spectrum.intensity.as_ref().iter())
                    .map(|(&mass, &intensity)| {
                        let mass = (mass - PROTON) * 1.0;
                        Peak {
                            mass,
                            intensity: intensity.to_f32(),
                        }
                    })
                    .collect::<Vec<_>>();
                (peaks, vec![])
            }
        };

        // process_ms2 deisotope branch already sorts; sort (or verify) the
        // other paths here -- covers MS1 spectra and MS2 with deisotope=false,
        // both of which return an empty peak_charges (see process_ms2/the `_`
        // arm above).
        if peak_charges.is_empty() {
            if self.assume_sorted_peaks {
                assert!(
                    is_sorted_by_mass(&peaks),
                    "Scan {}: `assume_sorted_peaks=true` but peaks are not \
                     mass-sorted ascending -- the caller's guarantee that \
                     upstream spectra arrive pre-sorted was violated. \
                     Continuing would silently break \
                     `select_matched_peaks`'s binary search precondition \
                     and produce wrong or missing fragment matches downstream, \
                     so this refuses to continue instead. Set \
                     `assume_sorted_peaks: false` to sort in-process instead.",
                    spectrum.id,
                );
            } else {
                peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
            }
        }
        // If peak_charges is non-empty the pairs were already sorted inside process_ms2

        let total_ion_current = peaks.iter().map(|peak| peak.intensity).sum::<f32>();

        ProcessedSpectrum {
            level: spectrum.ms_level,
            id: spectrum.id,
            file_id: spectrum.file_id,
            scan_start_time: spectrum.scan_start_time,
            ion_injection_time: spectrum.ion_injection_time,
            precursors: spectrum.precursors,
            peaks: peaks.into(),
            peak_charges,
            total_ion_current,
        }
    }

    pub fn process_with_mobility(&self, spectrum: RawSpectrum) -> ProcessedSpectrum<Vec<IMPeak>> {
        assert!(
            spectrum.ms_level == 1,
            "Logic error, mobility processing should only be used for MS1"
        );
        let mut peaks = spectrum
            .mz
            .iter()
            .zip(
                spectrum
                    .intensity
                    .iter()
                    .zip(spectrum.mobility.unwrap().iter()),
            )
            .map(|(&mass, (&intensity, &mobility))| {
                let mass = (mass - PROTON) * 1.0;
                IMPeak {
                    mass,
                    intensity,
                    mobility,
                }
            })
            .collect::<Vec<_>>();

        peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
        let total_ion_current = peaks.iter().map(|peak| peak.intensity).sum::<f32>();

        ProcessedSpectrum {
            level: spectrum.ms_level,
            id: spectrum.id,
            file_id: spectrum.file_id,
            scan_start_time: spectrum.scan_start_time,
            ion_injection_time: spectrum.ion_injection_time,
            precursors: spectrum.precursors,
            peaks,
            peak_charges: vec![],
            total_ion_current,
        }
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn columnar_peak_selection_matches_reference_including_ties() {
        let mut peaks = vec![
            Peak {
                mass: 99.0,
                intensity: 1000.0,
            },
            Peak {
                mass: 100.0,
                intensity: 5.0,
            },
            Peak {
                mass: 100.0,
                intensity: 10.0,
            },
            Peak {
                mass: 101.0,
                intensity: 10.0,
            },
            Peak {
                mass: 102.0,
                intensity: 0.0,
            },
            Peak {
                mass: 103.0,
                intensity: f32::NAN,
            },
        ];
        peaks.sort_by(|a, b| a.mass.total_cmp(&b.mass));
        let columns: PeakColumns = peaks.clone().into();
        for center in [0.0, 99.0, 100.0, 101.0, 102.0, 103.0, 1000.0] {
            for tolerance in [
                Tolerance::Da(0.0, 0.0),
                Tolerance::Da(-1.0, 1.0),
                Tolerance::Ppm(-10.0, 10.0),
            ] {
                for offset in [None, Some(-PROTON)] {
                    let (lo, hi) = tolerance.bounds(center);
                    let (lo, hi) = (
                        lo + offset.unwrap_or_default(),
                        hi + offset.unwrap_or_default(),
                    );
                    let mut best_intense = None;
                    let mut max_int = 0.0;
                    let mut best_closest = None;
                    let mut min_dist = f32::INFINITY;
                    let mut closest_int = 0.0;
                    for (idx, peak) in peaks.iter().enumerate() {
                        if peak.mass < lo || peak.mass > hi {
                            continue;
                        }
                        if peak.intensity >= max_int {
                            max_int = peak.intensity;
                            best_intense = Some(idx);
                        }
                        let dist = (peak.mass - center).abs();
                        if dist < min_dist || (dist == min_dist && peak.intensity >= closest_int)
                        {
                            min_dist = dist;
                            closest_int = peak.intensity;
                            best_closest = Some(idx);
                        }
                    }
                    let expected = match (best_intense, best_closest) {
                        (Some(most_intense), Some(closest)) => Some(PeakMatch {
                            most_intense,
                            closest,
                        }),
                        _ => None,
                    };
                    assert_eq!(
                        select_matched_peaks(&columns, center, tolerance, offset),
                        expected
                    );
                }
            }
        }
        assert_eq!(
            select_matched_peaks(&columns, 100.5, Tolerance::Da(-0.5, 0.5), None),
            Some(PeakMatch {
                most_intense: 3,
                closest: 3
            })
        );
        // Window containing the intense-but-farther 99.0 peak and the
        // weaker-but-closer 100.0 peaks -- proves the two fields are
        // genuinely independent, not always the same peak.
        assert_eq!(
            select_matched_peaks(&columns, 99.6, Tolerance::Da(-0.7, 0.7), None),
            Some(PeakMatch {
                most_intense: 0,
                closest: 2
            })
        );
        assert_eq!(
            select_matched_peaks(
                &PeakColumns::default(),
                100.0,
                Tolerance::Da(-1.0, 1.0),
                None
            ),
            None
        );
    }

    #[test]
    fn retain_peaks_keeps_all_columns_aligned() {
        let mut spectrum = ProcessedSpectrum {
            peaks: vec![
                Peak {
                    mass: 100.0,
                    intensity: 10.0,
                },
                Peak {
                    mass: 101.0,
                    intensity: 20.0,
                },
                Peak {
                    mass: 101.0,
                    intensity: 20.0,
                },
                Peak {
                    mass: 103.0,
                    intensity: 40.0,
                },
            ]
            .into(),
            peak_charges: vec![3, 1, 2, 2],
            ..Default::default()
        };
        spectrum.retain_peaks(|p| p.mass != 101.0);
        assert_eq!(spectrum.peaks.masses(), &[100.0, 103.0]);
        assert_eq!(spectrum.peaks.intensities(), &[10.0, 40.0]);
        assert_eq!(spectrum.peak_charges, vec![3, 2]);
        spectrum.retain_peaks(|_| false);
        assert_eq!(spectrum.peaks.len(), 0);
        assert!(spectrum.peak_charges.is_empty());
    }

    #[test]
    fn effective_precursor_tol_uses_own_window_when_set() {
        let p = Precursor {
            isolation_window: Some(Tolerance::Ppm(-5.0, 5.0)),
            ..Default::default()
        };
        let default = Tolerance::Ppm(-20.0, 20.0);
        match p.effective_precursor_tol(default) {
            Tolerance::Ppm(lo, hi) => {
                assert_eq!(lo, -5.0);
                assert_eq!(hi, 5.0);
            }
            other => panic!("expected Ppm, got {other:?}"),
        }
    }

    #[test]
    fn effective_precursor_tol_falls_back_to_default_when_unset() {
        let p = Precursor {
            isolation_window: None,
            ..Default::default()
        };
        let default = Tolerance::Ppm(-20.0, 20.0);
        match p.effective_precursor_tol(default) {
            Tolerance::Ppm(lo, hi) => {
                assert_eq!(lo, -20.0);
                assert_eq!(hi, 20.0);
            }
            other => panic!("expected Ppm, got {other:?}"),
        }
    }

    #[test]
    fn test_deisotope() {
        let mut mz = [
            800.9,
            800.9 + NEUTRON * 1.0,
            800.9 + NEUTRON * 2.0,
            803.4080,
            804.4108,
            805.4106,
            806.4116,
            810.0,
            812.0,
            812.0 + NEUTRON / 2.0,
        ];
        let mut int = [2., 1.5, 1., 4., 3., 2., 1., 1., 9.0, 4.5];
        let mut peaks = deisotope(&mut mz, &mut int, 2, 5.0, 800.91);

        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 2.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 1.0,
                    intensity: 2.5,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 2.0,
                    intensity: 1.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 6.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 3.0,
                    charge: Some(1),
                    envelope: Some(4),
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 1.0,
                    charge: Some(1),
                    envelope: Some(5),
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 4.5,
                    charge: Some(2),
                    envelope: Some(8),
                }
            ]
        );

        path_compression(&mut peaks);
        assert_eq!(
            peaks,
            vec![
                Deisotoped {
                    mz: 800.9,
                    intensity: 2.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 1.0,
                    intensity: 2.5,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 800.9 + NEUTRON * 2.0,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(1),
                },
                Deisotoped {
                    mz: 803.4080,
                    intensity: 10.0,
                    charge: Some(1),
                    envelope: None,
                },
                Deisotoped {
                    mz: 804.4108,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 805.4106,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 806.4116,
                    intensity: 0.0,
                    charge: Some(1),
                    envelope: Some(3),
                },
                Deisotoped {
                    mz: 810.0,
                    intensity: 1.0,
                    charge: None,
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0,
                    intensity: 13.5,
                    charge: Some(2),
                    envelope: None,
                },
                Deisotoped {
                    mz: 812.0 + NEUTRON / 2.0,
                    intensity: 0.0,
                    charge: Some(2),
                    envelope: Some(8),
                }
            ]
        );
    }

    fn mk_raw_spectrum(ms_level: u8, mz: Vec<f32>) -> RawSpectrum {
        let intensity = vec![1.0; mz.len()];
        RawSpectrum {
            ms_level,
            representation: Representation::Centroid,
            intensity,
            mz,
            ..RawSpectrum::default_with_file_id(0)
        }
    }

    #[test]
    fn borrowed_u32_processing_matches_owned_f32_including_rounding_ties() {
        let mz = [
            100.0,
            500.0,
            500.0 + NEUTRON,
            600.0,
            600.0 + NEUTRON / 2.0,
            700.0,
        ];
        let intensity = [0u32, 16_777_217, 16_777_216, 100, 25, u32::MAX];
        for count in [0, 1, mz.len()] {
            for deisotope in [false, true] {
                for top_n in [0, 1, 3, usize::MAX] {
                    let processor = SpectrumProcessor::new(top_n, deisotope, 0.0, false);
                    let precursors = vec![Precursor {
                        charge: Some(2),
                        ..Default::default()
                    }];
                    let owned: RawSpectrum = RawSpectrum {
                        mz: mz[..count].to_vec(),
                        intensity: intensity[..count].iter().map(|&x| x as f32).collect(),
                        ms_level: 2,
                        representation: Representation::Centroid,
                        precursors: precursors.clone(),
                        ..Default::default()
                    };
                    let borrowed = RawSpectrum {
                        mz: &mz[..count],
                        intensity: &intensity[..count],
                        ms_level: 2,
                        representation: Representation::Centroid,
                        precursors,
                        ..Default::default()
                    };
                    let expected = processor.process(owned);
                    let actual = processor.process(borrowed);
                    assert_eq!(actual.peaks, expected.peaks);
                    assert_eq!(actual.peak_charges, expected.peak_charges);
                    assert_eq!(
                        actual.total_ion_current.to_bits(),
                        expected.total_ion_current.to_bits()
                    );
                    assert!(actual.peaks.len() <= count.min(top_n));
                }
            }
        }
        let peaks = deisotope(&mz, &intensity, 2, 10.0, 0.0);
        assert_eq!(
            peaks[2].envelope, None,
            "rounded-equal intensities must not merge"
        );
        assert_eq!(
            peaks[4].envelope,
            Some(3),
            "ordinary isotope envelope must merge"
        );
    }

    #[test]
    fn is_sorted_by_mass_sorted_slice_is_true() {
        let peaks = vec![
            Peak {
                mass: 1.0,
                intensity: 1.0,
            },
            Peak {
                mass: 2.0,
                intensity: 1.0,
            },
            Peak {
                mass: 2.0,
                intensity: 1.0,
            }, // ties are fine (<=)
            Peak {
                mass: 3.0,
                intensity: 1.0,
            },
        ];
        assert!(is_sorted_by_mass(&peaks));
    }

    #[test]
    fn is_sorted_by_mass_unsorted_slice_is_false() {
        let peaks = vec![
            Peak {
                mass: 2.0,
                intensity: 1.0,
            },
            Peak {
                mass: 1.0,
                intensity: 1.0,
            },
        ];
        assert!(!is_sorted_by_mass(&peaks));
    }

    #[test]
    fn is_sorted_by_mass_empty_and_singleton_are_true() {
        assert!(is_sorted_by_mass(&[]));
        assert!(is_sorted_by_mass(&[Peak {
            mass: 5.0,
            intensity: 1.0
        }]));
    }

    #[test]
    fn process_deisotope_false_sorts_by_default() {
        let raw = mk_raw_spectrum(2, vec![500.0, 300.0, 400.0]);
        let sp = SpectrumProcessor::new(100, false, 0.0, false);
        let processed = sp.process(raw);
        assert!(is_sorted_by_mass(
            &processed.peaks.iter().collect::<Vec<_>>()
        ));
        assert_eq!(processed.peaks.len(), 3);
    }

    #[test]
    fn process_assume_sorted_peaks_true_passes_through_sorted_input() {
        let raw = mk_raw_spectrum(2, vec![300.0, 400.0, 500.0]);
        let sp = SpectrumProcessor::new(100, false, 0.0, true);
        let processed = sp.process(raw);
        assert!(is_sorted_by_mass(
            &processed.peaks.iter().collect::<Vec<_>>()
        ));
        assert_eq!(processed.peaks.len(), 3);
    }

    #[test]
    #[should_panic(expected = "not mass-sorted ascending")]
    fn process_assume_sorted_peaks_true_panics_on_unsorted_ms2_input() {
        let raw = mk_raw_spectrum(2, vec![500.0, 300.0, 400.0]);
        let sp = SpectrumProcessor::new(100, false, 0.0, true);
        sp.process(raw);
    }

    #[test]
    #[should_panic(expected = "not mass-sorted ascending")]
    fn process_assume_sorted_peaks_true_panics_on_unsorted_ms1_input() {
        let raw = mk_raw_spectrum(1, vec![500.0, 300.0, 400.0]);
        let sp = SpectrumProcessor::new(100, false, 0.0, true);
        sp.process(raw);
    }

    #[test]
    fn process_assume_sorted_peaks_true_is_noop_when_deisotope_true() {
        // deisotope=true's branch never reaches the assume_sorted_peaks
        // check -- its own sort runs first, for a different purpose
        // (selecting top-N by intensity) -- so this must not panic despite
        // deliberately unsorted input.
        let raw = mk_raw_spectrum(2, vec![500.0, 300.0, 400.0]);
        let sp = SpectrumProcessor::new(100, true, 0.0, true);
        let processed = sp.process(raw);
        assert!(is_sorted_by_mass(
            &processed.peaks.iter().collect::<Vec<_>>()
        ));
    }
}
