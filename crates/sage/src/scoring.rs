use crate::database::{IndexedDatabase, IndexedQuery, PeptideIx};
use crate::heap::bounded_min_heapify;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::{Tolerance, NEUTRON, PROTON};
use crate::ms2_similarity::{self, N_FRAGMENT_SLOTS};
use crate::spectrum::{Peak, Precursor, ProcessedSpectrum};
use half::f16;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::ops::AddAssign;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};

thread_local! {
    static WINDOWS_SCRATCH: RefCell<Vec<(f32, Tolerance)>> = RefCell::new(Vec::new());
}

#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum ScoreType {
    SageHyperScore,
    OpenMSHyperScore,
}

/// Which key `build_features` ranks/truncates candidates on. Runtime-selectable
/// (config `ranking_score`, same `Option<T>` + `Input::build()`-default shape as
/// `score_type` above) rather than only ever compiled in — `combined_score` is
/// always computed regardless of this setting (cheap: dense-array z² lookups),
/// just not always used for ranking.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub enum RankingScore {
    /// Raw X!Tandem hyperscore — SAGE's ranking behavior before external
    /// RT/IIM predictions could influence candidate ranking at all.
    Hyperscore,
    /// `hyperscore - 0.5 * (z_rt_external² + z_iim_external²)` — see `Score::combined_score`.
    CombinedScore,
}

/// Structure to hold temporary scores
#[derive(Copy, Clone, Default, Debug, PartialEq, PartialOrd)]
struct Score {
    peptide: PeptideIx,
    matched_b: u16,
    matched_y: u16,
    summed_b: f32,
    summed_y: f32,
    longest_b: usize,
    longest_y: usize,
    hyperscore: f64,
    /// `hyperscore - 0.5 * (z_rt_external² + z_iim_external²)` — the actual
    /// ranking/retention key (`build_features` sorts and truncates to
    /// `report_psms` on this, not raw `hyperscore`). Equal to `hyperscore`
    /// whenever `--predicted-rt`/`--predicted-iim` aren't configured (the
    /// z² terms are 0.0 in that case) — see `external_z2`.
    combined_score: f64,
    ppm_difference: f32,
    precursor_charge: u8,
    isotope_error: i8,
    /// MS2 predicted-vs-observed fragment intensity similarity metrics
    /// (`ms2_similarity` module, MSBooster parity -- see
    /// `docs/ai/predicted_fragment_intensity.md`) -- reported `Feature`
    /// fields only, never part of `combined_score`/ranking (contrast the
    /// RT/IIM z² terms). All `0.0` (or `-1.0` for the two correlations,
    /// their own sentinel) when `--predicted-fragment-intensity-*` isn't
    /// configured, or this candidate's `(peptide, charge)` has no cache
    /// entry.
    ms2_entropy_similarity: f32,
    ms2_weighted_entropy_similarity: f32,
    ms2_heuristic_entropy_similarity: f32,
    ms2_cosine_similarity: f32,
    ms2_dot_product: f32,
    ms2_spectral_contrast_angle: f32,
    ms2_euclidean_similarity: f32,
    ms2_bray_curtis_similarity: f32,
    ms2_pearson_corr: f32,
    ms2_spearman_corr: f32,
    ms2_hypergeometric_probability: f32,
    ms2_intersection: u32,
    ms2_top6_matched_intensity: f32,
}

impl Eq for Score {}

impl Ord for Score {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.hyperscore
            .partial_cmp(&other.hyperscore)
            .unwrap_or(std::cmp::Ordering::Less)
    }
}

/// Preliminary score - # of matched peaks for each candidate peptide
#[derive(Copy, Clone, Default, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct PreScore {
    matched: u16,
    peptide: PeptideIx,
    precursor_charge: u8,
    isotope_error: i8,
}

/// Store preliminary scores & stats for first pass search for a query spectrum
#[derive(Clone, Default)]
struct InitialHits {
    matched_peaks: usize,
    // Number of peptide candidates with > 0 matched peaks
    scored_candidates: usize,
    preliminary: Vec<PreScore>,
}

impl AddAssign<InitialHits> for InitialHits {
    fn add_assign(&mut self, rhs: InitialHits) {
        self.matched_peaks += rhs.matched_peaks;
        self.scored_candidates += rhs.scored_candidates;

        self.preliminary.extend(rhs.preliminary);
    }
}

#[derive(Serialize, Clone, Debug, Default)]
/// Features of a candidate peptide spectrum match
pub struct Feature {
    #[serde(skip_serializing)]
    pub peptide_idx: PeptideIx,
    // psm_id help to match with matched fragments table.
    pub psm_id: usize,
    pub peptide_len: usize,
    /// Spectrum id
    pub spec_id: String,
    /// File identifier
    pub file_id: usize,
    /// PSM rank
    pub rank: u32,
    /// Target/Decoy label, -1 is decoy, 1 is target
    pub label: i32,
    /// Experimental mass
    pub expmass: f32,
    /// Calculated mass
    pub calcmass: f32,
    /// Reported precursor charge
    pub charge: u8,
    /// Retention time
    pub rt: f32,
    /// Globally aligned retention time
    pub aligned_rt: f32,
    /// Predicted RT, if enabled
    pub predicted_rt: f32,
    /// Difference between predicted & observed RT
    pub delta_rt_model: f32,
    /// Externally-predicted RT (`--predicted-rt`, e.g. Chronologer,
    /// calibrated per-run), independent of `predicted_rt`/`delta_rt_model`
    /// above (SAGE's own in-run composition-regression model — kept as a
    /// separate feature, not replaced). 0.0 when `--predicted-rt` isn't
    /// configured. See `plans/lda_external_rt_iim_features.md`.
    pub predicted_rt_external: f32,
    /// `((aligned_rt - predicted_rt_external) / rt_sigma)²` — a z² "badness"
    /// feature for SAGE's LDA (`ml/linear_discriminant.rs`). 0.0 when
    /// `--predicted-rt` isn't configured (and then not included in the LDA
    /// at all — see `plans/lda_external_rt_iim_features.md`).
    pub delta_rt_z2_external: f32,
    /// Ion mobility
    pub ims: f32,
    /// Predicted ion mobility, if enabled
    pub predicted_ims: f32,
    /// Difference between predicted & observed ion mobility
    pub delta_ims_model: f32,
    /// Externally-predicted IIM (`--predicted-iim`, e.g. IM2Deep),
    /// independent of `predicted_ims`/`delta_ims_model` above (SAGE's own
    /// in-run model). 0.0 when `--predicted-iim` isn't configured.
    pub predicted_ims_external: f32,
    /// `((ims - predicted_ims_external) / iim_sigma)²`, same shape as
    /// `delta_rt_z2_external`. 0.0 when `--predicted-iim` isn't configured.
    pub delta_ims_z2_external: f32,
    /// Difference between expmass and calcmass
    pub delta_mass: f32,
    /// C13 isotope error
    pub isotope_error: f32,
    /// Average ppm delta mass for matched fragments
    pub average_ppm: f32,
    /// X!Tandem hyperscore
    pub hyperscore: f64,
    /// Difference between hyperscore of this candidate, and the next best candidate
    pub delta_next: f64,
    /// Difference between hyperscore of this candidate, and the best candidate
    pub delta_best: f64,
    /// Number of matched theoretical fragment ions
    pub matched_peaks: u32,
    /// Longest b-ion series
    pub longest_b: u32,
    /// Longest y-ion series
    pub longest_y: u32,
    /// Longest y-ion series, divided by peptide length
    pub longest_y_pct: f32,
    /// Number of missed cleavages
    pub missed_cleavages: u8,
    /// Fraction of matched MS2 intensity
    pub matched_intensity_pct: f32,
    /// Number of scored candidates for this spectrum
    pub scored_candidates: u32,
    /// Probability of matching exactly N peaks across all candidates Pr(x=k)
    pub poisson: f64,
    /// Combined score from linear discriminant analysis, used for FDR calc
    pub discriminant_score: f32,
    /// Posterior error probability for this PSM / local FDR
    pub posterior_error: f32,
    /// Assigned q_value
    pub spectrum_q: f32,
    pub peptide_q: f32,
    pub protein_q: f32,
    pub protein_group_q: f32,

    pub ms2_intensity: f32,

    /// MS2 predicted-vs-observed fragment intensity similarity metrics
    /// (MSBooster parity, see `ms2_similarity` module and
    /// `docs/ai/predicted_fragment_intensity.md`) -- all `0.0` (or `-1.0`
    /// for the two correlations, their own sentinel) when
    /// `--predicted-fragment-intensity-*` isn't configured, or this
    /// `(peptide, charge)` has no entry in the cache. Feature-only, no hard
    /// filtering (contrast `predicted_rt_external`/`delta_rt_z2_external`,
    /// which also feed `evict_rt_iim_mismatches`/`combined_score`).
    pub ms2_entropy_similarity: f32,
    pub ms2_weighted_entropy_similarity: f32,
    pub ms2_heuristic_entropy_similarity: f32,
    pub ms2_cosine_similarity: f32,
    pub ms2_dot_product: f32,
    pub ms2_spectral_contrast_angle: f32,
    pub ms2_euclidean_similarity: f32,
    pub ms2_bray_curtis_similarity: f32,
    pub ms2_pearson_corr: f32,
    pub ms2_spearman_corr: f32,
    pub ms2_hypergeometric_probability: f32,
    pub ms2_intersection: u32,
    pub ms2_top6_matched_intensity: f32,

    pub protein_groups: Option<String>,
    pub num_protein_groups: u32,

    pub fragments: Option<Fragments>,
}

/// Matching Fragment details
#[derive(Serialize, Default, Clone, Debug)]
pub struct Fragments {
    #[serde(skip_serializing)]
    pub charges: Vec<i32>,
    pub kinds: Vec<Kind>,
    pub fragment_ordinals: Vec<i32>,
    pub intensities: Vec<f32>,
    pub mz_calculated: Vec<f32>,
    pub mz_experimental: Vec<f32>,
}

static PSM_COUNTER: AtomicUsize = AtomicUsize::new(1);

fn increment_psm_counter() -> usize {
    PSM_COUNTER.fetch_add(1, Ordering::Relaxed)
}

/// Stirling's approximation for log factorial
fn lnfact(n: u16) -> f64 {
    if n == 0 {
        1.0
    } else {
        let n = n as f64;
        n * n.ln() - n + 0.5 * n.ln() + 0.5 * (std::f64::consts::PI * 2.0 * n).ln()
    }
}

impl ScoreType {
    pub fn score(&self, matched_b: u16, matched_y: u16, summed_b: f32, summed_y: f32) -> f64 {
        let score = match self {
            // Calculate the X!Tandem hyperscore
            Self::SageHyperScore => {
                let i = (summed_b + 1.0) as f64 * (summed_y + 1.0) as f64;

                i.ln() + lnfact(matched_b) + lnfact(matched_y)
            }
            // Calculate the OpenMS flavour hyperscore
            Self::OpenMSHyperScore => {
                let summed_intensity = summed_b + summed_y;

                summed_intensity.ln_1p() as f64 + lnfact(matched_b) + lnfact(matched_y)
            }
        };
        if score.is_finite() {
            score
        } else {
            255.0
        }
    }
}

impl Score {
    /// Calculate the hyperscore for a given PSM choosing between implementations based on `score_type`
    fn hyperscore(&self, score_type: ScoreType) -> f64 {
        score_type.score(self.matched_b, self.matched_y, self.summed_b, self.summed_y)
    }

    /// The ranking/retention key `build_features` actually sorts and
    /// truncates candidates on, per `Scorer::ranking_score`.
    fn rank_key(&self, ranking_score: RankingScore) -> f64 {
        match ranking_score {
            RankingScore::Hyperscore => self.hyperscore,
            RankingScore::CombinedScore => self.combined_score,
        }
    }
}

pub struct Scorer<'db> {
    pub db: &'db IndexedDatabase,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    /// What is the minimum number of matched b and y ion peaks to report PSMs for?
    pub min_matched_peaks: u16,
    /// Precursor isotope error lower bounds (e.g. -1)
    pub min_isotope_err: i8,
    /// Precursor isotope error upper bounds (e.g. 3)
    pub max_isotope_err: i8,
    pub min_precursor_charge: u8,
    pub max_precursor_charge: u8,
    pub override_precursor_charge: bool,
    pub max_fragment_charge: Option<u8>,
    pub chimera: bool,
    pub report_psms: usize,

    // Rather than use a fixed precursor tolerance, dynamically alter
    // the precursor tolerance window based on MS2 isolation window and charge
    pub wide_window: bool,
    pub annotate_matches: bool,
    pub score_type: ScoreType,
    pub ranking_score: RankingScore,

    /// Externally-predicted RT, indexed by peptide index (`db.peptides`),
    /// used to evict candidates whose predicted RT falls outside `rt_tol`
    /// of the observed spectrum's RT. `rt_tol` is required (validated at
    /// config-load time) whenever this is `Some`. Independent of
    /// `predicted_iim` — either, both, or neither may be set. Resolved
    /// once from a `sequence -> rt` map by `Runner::resolve_predicted_rt`
    /// (not read directly by `Scorer`) — the index form avoids calling
    /// `Peptide::to_string()` per (candidate, spectrum) during search, a
    /// real measured ~37% `run_sage` overhead before this. See
    /// `plans/rt_iim_independent_dimensions.md`.
    pub predicted_rt: Option<&'db [Option<f32>]>,
    /// Externally-predicted IIM, dense-indexed via [`iim_dense_slot`] (one
    /// slot per `(peptide index, charge)` combination in
    /// `[min_precursor_charge, max_precursor_charge]`), used to evict
    /// candidates whose predicted IIM falls outside `mobility_tol` of the
    /// observed spectrum's IIM. `mobility_tol` is required (validated at
    /// config-load time) whenever this is `Some`. Independent of
    /// `predicted_rt`. Resolved once by `Runner::resolve_predicted_iim`,
    /// same reasoning as `predicted_rt`. A dense array beats a
    /// `HashMap<(usize, u8), f32>` here: a standalone benchmark at real
    /// F9477 scale (18.9M entries, 20M random-access queries, fixed seed)
    /// measured 239ns/query for `std::HashMap` vs 10.4ns/query for a dense
    /// array — 23x, and smaller too (no per-entry key storage or
    /// hashmap load-factor overhead). See
    /// `plans/rt_iim_independent_dimensions.md`.
    pub predicted_iim: Option<&'db [Option<f32>]>,
    /// Robust (MAD-based) scale of the `predicted_rt` residual as a
    /// function of observed `scan_start_time`, already converted to minutes
    /// (`sage-cli/src/input.rs`'s `spline_secs_to_minutes` at config-load
    /// time). Normalizes `Feature::delta_rt_z2_external` into a z² LDA
    /// feature, evaluated at each candidate's own `scan_start_time` rather
    /// than a single global scale — see
    /// `plans/lda_external_rt_iim_features.md` and
    /// `plans/rt_heteroscedastic_tolerance_spline.md`.
    pub rt_sigma: Option<crate::spline::LinearSpline>,
    /// Same as `rt_sigma`, for `predicted_iim` — unitless (1/K0), no
    /// conversion needed.
    pub iim_sigma: Option<f32>,
    /// RT tolerance window as a function of observed `scan_start_time`,
    /// already in minutes (converted from the config's `rt_tol_sec` at load
    /// time) to match `ProcessedSpectrum::scan_start_time`'s own unit. A
    /// flat window is just a 2-node spline with identical values at both
    /// nodes — there is no separate flat-tolerance representation.
    pub rt_tol: Option<crate::spline::ValueTolSpline>,
    /// IIM tolerance window as a function of observed
    /// `Precursor::inverse_ion_mobility`.
    pub mobility_tol: Option<crate::spline::ValueTolSpline>,

    /// Externally-predicted MS2 fragment intensities
    /// (`--predicted-fragment-intensity-index`/`-cache`, Prosit `compact_trt`
    /// via `git/featureprediction`) -- **feature-only, no hard filtering**
    /// (contrast `predicted_rt`/`predicted_iim` above): used solely to
    /// compute `Feature::ms2_entropy_similarity`, never to evict or re-rank
    /// candidates. Dense-indexed by `(peptide index, charge)` via the same
    /// [`iim_dense_slot`] shape `predicted_iim` uses -- each entry is a
    /// `(start, end)` half-open row range into
    /// `predicted_fragment_intensity_annotation_id`/`_intensity` below.
    /// Resolved once by `Runner::resolve_predicted_fragment_intensity`, same
    /// reasoning as `predicted_rt`/`predicted_iim` (avoids `Peptide::
    /// to_string()` per candidate). See
    /// `docs/ai/predicted_fragment_intensity.md`.
    pub predicted_fragment_intensity_index: Option<&'db [Option<(u64, u64)>]>,
    /// The shared, ever-growing `arrays.mmappet` annotation-id column
    /// (`git/featureprediction`'s SAGE-order fragment-slot numbering, see
    /// `ms2_similarity::fragment_annotation_id`), sliced per candidate via
    /// `predicted_fragment_intensity_index`'s `(start, end)`.
    pub predicted_fragment_intensity_annotation_id: Option<&'db [u8]>,
    /// The shared `arrays.mmappet` predicted-intensity column (native FP16
    /// from Prosit `compact_trt`), parallel to
    /// `predicted_fragment_intensity_annotation_id`.
    pub predicted_fragment_intensity: Option<&'db [f16]>,
}

#[inline(always)]
/// Calculate upper bound (excluded) of the charge state range to use for
/// searching fragment ions (1..N)
/// If user has configured max_fragment_charge, potentially override precursor
/// charge
fn max_fragment_charge(max_fragment_charge: Option<u8>, precursor_charge: u8) -> u8 {
    precursor_charge
        .min(
            max_fragment_charge
                .map(|c| c + 1)
                .unwrap_or(precursor_charge),
        )
        .max(2)
}

impl<'db> Scorer<'db> {
    /// Perform a quick first-pass scoring, where we consider a peptide "identified"
    /// if it meets the following criterion:
    ///  * prefilter_low_memory = true: in the top `report_psms` hits for a spectrum
    ///  * prefilter_low_memory = false: has at least `min_matched_peaks` fragment ion matches
    /// * `keep`: A vector of atomic bools is used to maintain an identification list across scans
    pub fn quick_score(
        &self,
        query: &ProcessedSpectrum<Peak>,
        prefilter_low_memory: bool,
        keep: &[AtomicBool],
    ) {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });
        let hits = self.initial_hits(query, precursor);

        if prefilter_low_memory {
            let mut score_vector = hits
                .preliminary
                .iter()
                .filter_map(|pre| {
                    if pre.peptide == PeptideIx::default() {
                        return None;
                    }
                    let (score, _) = self.score_candidate(query, pre);
                    if (score.matched_b + score.matched_y) < self.min_matched_peaks {
                        return None;
                    }
                    Some(score)
                })
                .collect::<Vec<_>>();

            let k = self.report_psms.min(score_vector.len());
            bounded_min_heapify(&mut score_vector, k);
            for score in &score_vector[..k] {
                keep[score.peptide.0 as usize].store(true, Ordering::Relaxed);
            }
        } else {
            for pre in &hits.preliminary {
                if pre.peptide != PeptideIx::default() {
                    keep[pre.peptide.0 as usize].store(true, Ordering::Relaxed);
                }
            }
        }
    }

    pub fn score(&self, query: &ProcessedSpectrum<crate::spectrum::Peak>) -> Vec<Feature> {
        assert_eq!(
            query.level, 2,
            "internal bug, trying to score a non-MS2 scan!"
        );
        match self.chimera {
            true => self.score_chimera_fast(query),
            false => self.score_standard(query),
        }
    }

    /// Perform a k-select and truncation of an [`InitialHits`] list.
    ///
    /// Determine how many candidates to actually calculate hyperscore for.
    /// Hyperscore is relatively computationally expensive, so we don't want
    /// to calculate it for every possible candidate (100s - 10,000s depending on search)
    /// when we are only going to report a few PSMs. But we also want to calculate
    /// it for enough candidates that we don't accidentally miss the best hit!
    ///
    /// Given that hyperscore is dominated by the number of matched peaks, it seems
    /// reasonable to assume that the highest hyperscore will belong to one of the
    /// top 50 candidates sorted by # of matched peaks.
    fn trim_hits(&self, hits: &mut InitialHits) {
        let k = 50.clamp(
            (self.report_psms * 2).min(hits.preliminary.len()),
            hits.preliminary.len(),
        );
        bounded_min_heapify(&mut hits.preliminary, k);
        hits.preliminary.truncate(k);
    }

    /// Preliminary Score, return # of matched peaks per candidate
    /// Returned hits are guaranteed to be the top-K hits (see above comment)
    /// from among all potential candidates, but the returned vector is not
    /// in sorted order.
    fn matched_peaks_with_isotope(
        &self,
        query: &ProcessedSpectrum<crate::spectrum::Peak>,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
        isotope_error: i8,
    ) -> InitialHits {
        let candidates = self
            .db
            .query(precursor_mass - isotope_error as f32 * NEUTRON, precursor_tol);

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, precursor_charge);
        // Allocate space for all potential candidates - many potential candidates
        let potential = candidates.pre_idx_hi - candidates.pre_idx_lo + 1;
        let mut hits = InitialHits {
            matched_peaks: 0,
            scored_candidates: 0,
            preliminary: vec![PreScore::default(); potential],
        };

        // Collect every (mass, fragment_tol) window for this whole
        // spectrum (every peak x every fragment charge) up front, then
        // search them in one batched call -- `page_search_batch` computes
        // each page's precursor-scoped inner range once and reuses it
        // across every window landing on that page, instead of redoing
        // that binary search independently per (peak, charge) the way a
        // `page_search` call per window did. `windows` is thread-local
        // scratch, reused across calls, not freshly allocated each time
        // (same reasoning as `page_search_batch`'s own internal scratch).
        // See `docs/ai/reuse_index_bins.md`.
        WINDOWS_SCRATCH.with(|windows| {
            let mut windows = windows.borrow_mut();
            windows.clear();
            for peak in query.peaks.iter() {
                for charge in 1..max_fragment_charge {
                    windows.push((peak.mass * charge as f32, self.fragment_tol));
                }
            }

            candidates.page_search_batch(&windows, |frag| {
                let idx = frag.peptide_index.0 as usize - candidates.pre_idx_lo;
                let sc = &mut hits.preliminary[idx];
                if sc.matched == 0 {
                    hits.scored_candidates += 1;
                    sc.precursor_charge = precursor_charge;
                    sc.peptide = frag.peptide_index;
                    sc.isotope_error = isotope_error;
                }

                sc.matched += 1;
                hits.matched_peaks += 1;
            });
        });
        if hits.matched_peaks == 0 {
            return hits;
        }

        if self.predicted_rt.is_some() || self.predicted_iim.is_some() {
            self.evict_rt_iim_mismatches(query, &candidates, &mut hits);
        }

        self.trim_hits(&mut hits);
        hits
    }

    /// Maps `(peptide_idx, charge)` to a slot in a dense
    /// `predicted_iim`-shaped array of length `n_peptides *
    /// (max_charge - min_charge + 1)`, or `None` if `charge` falls outside
    /// `[min_charge, max_charge]` entirely (can't have an entry — same
    /// permissive treatment as a missing key would get from a map). Shared
    /// between `Runner::resolve_predicted_iim` (which builds the array)
    /// and `evict_rt_iim_mismatches` (which reads it) so the indexing
    /// arithmetic can't drift between the two call sites. See
    /// `plans/rt_iim_independent_dimensions.md`.
    pub fn iim_dense_slot(
        peptide_idx: usize,
        charge: u8,
        min_charge: u8,
        max_charge: u8,
    ) -> Option<usize> {
        if charge < min_charge || charge > max_charge {
            return None;
        }
        let charge_span = (max_charge - min_charge + 1) as usize;
        Some(peptide_idx * charge_span + (charge - min_charge) as usize)
    }

    /// Evict candidates (in place) whose predicted RT and/or IIM falls
    /// outside `rt_tol`/`mobility_tol` of `query`'s observed values. Each
    /// dimension is checked independently and only if its map
    /// (`self.predicted_rt`/`self.predicted_iim`) is configured at all — an
    /// unconfigured dimension never evicts anything. Candidates with no
    /// entry in a configured map for their `sequence`/`(sequence, charge)`
    /// are left alone on that dimension (permissive — avoids rejecting due
    /// to prediction-coverage gaps). Missing observed
    /// `inverse_ion_mobility` (non-PASEF data) skips the IIM check only;
    /// the RT check still applies.
    fn evict_rt_iim_mismatches(
        &self,
        query: &ProcessedSpectrum<crate::spectrum::Peak>,
        candidates: &IndexedQuery,
        hits: &mut InitialHits,
    ) {
        let rt_bounds = self.predicted_rt.map(|_| {
            self.rt_tol
                .as_ref()
                .expect("validated at config-load time: rt_tol required with predicted_rt")
                .tolerance_at(query.scan_start_time)
                .bounds(query.scan_start_time)
        });
        let iim_bounds = self.predicted_iim.and_then(|_| {
            query
                .precursors
                .first()
                .and_then(|p| p.inverse_ion_mobility)
                .map(|observed| {
                    self.mobility_tol
                        .as_ref()
                        .expect(
                            "validated at config-load time: mobility_tol required with predicted_iim",
                        )
                        .tolerance_at(observed)
                        .bounds(observed)
                })
        });

        for (i, sc) in hits.preliminary.iter_mut().enumerate() {
            if sc.matched == 0 {
                continue;
            }
            let peptide_idx = candidates.pre_idx_lo + i;

            let rt_ok = rt_bounds.map_or(true, |(lo, hi)| {
                self.predicted_rt
                    .and_then(|by_idx| by_idx[peptide_idx])
                    .map_or(true, |rt| rt >= lo && rt <= hi)
            });
            let iim_ok = iim_bounds.map_or(true, |(lo, hi)| {
                Self::iim_dense_slot(
                    peptide_idx,
                    sc.precursor_charge,
                    self.min_precursor_charge,
                    self.max_precursor_charge,
                )
                .and_then(|slot| self.predicted_iim.and_then(|dense| dense[slot]))
                .map_or(true, |iim| iim >= lo && iim <= hi)
            });

            if !rt_ok || !iim_ok {
                hits.matched_peaks -= sc.matched as usize;
                hits.scored_candidates -= 1;
                *sc = PreScore::default();
            }
        }
    }

    /// External (Chronologer/IM2Deep, per-run-calibrated) RT/IIM z² terms for
    /// one candidate — `((observed - predicted) / sigma)²`, 0.0 when that
    /// dimension isn't configured (`--predicted-rt`/`--predicted-iim` unset
    /// or `sigma <= 0.0`) or this peptide/charge has no entry in the dense
    /// array. Shared by `build_features`'s pre-sort `combined_score` pass and
    /// its post-sort `Feature`-building pass so the two can't drift apart.
    fn external_z2(
        &self,
        query: &ProcessedSpectrum<Peak>,
        peptide_idx: usize,
        charge: u8,
    ) -> (f32, f32, f32, f32) {
        // `rt_sigma` is now RT-dependent (a `LinearSpline`, not a bare
        // scalar), so it must be `.eval(scan_start_time)`-evaluated before
        // the `sigma > 0.0` guard can be checked -- can no longer gate on
        // that in the match pattern itself the way a `Copy` scalar could.
        let (predicted_rt_external, delta_rt_z2_external) =
            match (self.predicted_rt, self.rt_sigma.as_ref()) {
                (Some(by_idx), Some(sigma_spline)) => match by_idx[peptide_idx] {
                    Some(rt) => {
                        let sigma = sigma_spline.eval(query.scan_start_time);
                        if sigma > 0.0 {
                            let z = (query.scan_start_time - rt) / sigma;
                            (rt, z * z)
                        } else {
                            (0.0, 0.0)
                        }
                    }
                    None => (0.0, 0.0),
                },
                _ => (0.0, 0.0),
            };

        let observed_ims = query
            .precursors
            .first()
            .and_then(|p| p.inverse_ion_mobility);
        let (predicted_ims_external, delta_ims_z2_external) =
            match (self.predicted_iim, self.iim_sigma, observed_ims) {
                (Some(dense), Some(sigma), Some(observed)) if sigma > 0.0 => {
                    Self::iim_dense_slot(
                        peptide_idx,
                        charge,
                        self.min_precursor_charge,
                        self.max_precursor_charge,
                    )
                    .and_then(|slot| dense[slot])
                    .map_or((0.0, 0.0), |ims| {
                        let z = (observed - ims) / sigma;
                        (ims, z * z)
                    })
                }
                _ => (0.0, 0.0),
            };

        (
            predicted_rt_external,
            delta_rt_z2_external,
            predicted_ims_external,
            delta_ims_z2_external,
        )
    }

    fn matched_peaks(
        &self,
        query: &ProcessedSpectrum<Peak>,
        precursor_mass: f32,
        precursor_charge: u8,
        precursor_tol: Tolerance,
    ) -> InitialHits {
        if self.min_isotope_err != self.max_isotope_err {
            let mut hits = (self.min_isotope_err..=self.max_isotope_err).fold(
                InitialHits::default(),
                |mut hits, isotope| {
                    hits += self.matched_peaks_with_isotope(
                        query,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                        isotope,
                    );
                    hits
                },
            );
            self.trim_hits(&mut hits);
            hits
        } else {
            self.matched_peaks_with_isotope(
                query,
                precursor_mass,
                precursor_charge,
                precursor_tol,
                0,
            )
        }
    }

    fn initial_hits(&self, query: &ProcessedSpectrum<Peak>, precursor: &Precursor) -> InitialHits {
        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        // Search in wide-window/DIA mode
        let mut hits = if self.wide_window {
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    let precursor_tol = precursor
                        .isolation_window
                        .unwrap_or(Tolerance::Da(-2.4, 2.4))
                        * precursor_charge as f32;
                    hits +=
                        self.matched_peaks(query, precursor_mass, precursor_charge, precursor_tol);
                    hits
                },
            )
        } else if precursor.charge.is_some() && !self.override_precursor_charge {
            let charge = precursor.charge.unwrap();
            // Charge state is already annotated for this precusor, only search once
            let precursor_mass = mz * charge as f32;
            let precursor_tol = precursor.effective_precursor_tol(self.precursor_tol);
            self.matched_peaks(query, precursor_mass, charge, precursor_tol)
        } else {
            // Not all selected ion precursors have charge states annotated (or user has set
            // `override_precursor_charge`)
            // assume it could be z=2, z=3, z=4 and search all three
            let precursor_tol = precursor.effective_precursor_tol(self.precursor_tol);
            (self.min_precursor_charge..=self.max_precursor_charge).fold(
                InitialHits::default(),
                |mut hits, precursor_charge| {
                    let precursor_mass = mz * precursor_charge as f32;
                    hits += self.matched_peaks(
                        query,
                        precursor_mass,
                        precursor_charge,
                        precursor_tol,
                    );
                    hits
                },
            )
        };
        self.trim_hits(&mut hits);
        hits
    }

    /// Score a single [`ProcessedSpectrum`] against the database
    pub fn score_standard(&self, query: &ProcessedSpectrum<Peak>) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let hits = self.initial_hits(query, precursor);
        let mut features = Vec::with_capacity(self.report_psms);
        self.build_features(query, precursor, &hits, self.report_psms, &mut features);
        features
    }

    /// Given a set of [`InitialHits`] against a query spectrum, prepare N=`report_psms`
    /// best PSMs ([`Feature`])
    fn build_features(
        &self,
        query: &ProcessedSpectrum<Peak>,
        precursor: &Precursor,
        hits: &InitialHits,
        report_psms: usize,
        features: &mut Vec<Feature>,
    ) {
        let mut score_vector = hits
            .preliminary
            .iter()
            .filter(|score| score.peptide != PeptideIx::default())
            .map(|pre| self.score_candidate(query, pre))
            .filter(|s| (s.0.matched_b + s.0.matched_y) >= self.min_matched_peaks)
            .map(|(mut score, fragments)| {
                let (_, z_rt2, _, z_iim2) =
                    self.external_z2(query, score.peptide.0 as usize, score.precursor_charge);
                score.combined_score =
                    score.hyperscore - 0.5 * (z_rt2 as f64 + z_iim2 as f64);
                (score, fragments)
            })
            .collect::<Vec<_>>();

        // Ranking/retention key -- `Scorer::ranking_score` selects between
        // raw hyperscore and `combined_score` (hyperscore softly penalized by
        // external RT/IIM implausibility; 0 penalty, i.e. identical to plain
        // hyperscore, when `--predicted-rt`/`--predicted-iim` aren't
        // configured). Under `CombinedScore`, candidates closer to their
        // predicted RT/IIM rank higher among otherwise similar hyperscores,
        // rather than being hard-evicted.
        score_vector.sort_by(|a, b| {
            b.0.rank_key(self.ranking_score)
                .total_cmp(&a.0.rank_key(self.ranking_score))
        });

        // Expected value for poisson distribution
        // (average # of matches peaks/peptide candidate)
        let lambda = hits.matched_peaks as f64 / hits.scored_candidates as f64;

        // Sage operates on masses without protons; [M] instead of [MH+]
        let mz = precursor.mz - PROTON;

        for idx in 0..report_psms.min(score_vector.len()) {
            let score = score_vector[idx].0;
            let fragments: Option<Fragments> = score_vector[idx].1.take();
            let psm_id = increment_psm_counter();

            let peptide = &self.db[score.peptide];
            let precursor_mass = mz * score.precursor_charge as f32;

            // Margins are on `rank_key` (the actual ranking key per
            // `ranking_score`), not unconditionally raw `hyperscore` —
            // consistent with the sort above.
            let next = score_vector
                .get(idx + 1)
                .map(|score| score.0.rank_key(self.ranking_score))
                .unwrap_or_default();

            let best = score_vector
                .first()
                .map(|score| score.0.rank_key(self.ranking_score))
                .expect("we know that index 0 is valid");

            // Poisson distribution log10 probability mass function
            // Computed directly in log space to avoid overflow from lambda.powi(k)
            // log10(PMF) = (k*ln(lambda) - lambda - lnfact(k)) / ln(10)
            let k = score.matched_b + score.matched_y;
            let log10_poisson =
                (k as f64 * lambda.ln() - lambda - lnfact(k)) / std::f64::consts::LN_10;

            let isotope_error = score.isotope_error as f32 * NEUTRON;
            let delta_mass = (precursor_mass - peptide.monoisotopic - isotope_error) * 2E6
                / (precursor_mass - isotope_error + peptide.monoisotopic);

            // External (Chronologer/IM2Deep, per-run-calibrated) RT/IIM,
            // independent of SAGE's own in-run composition-regression model
            // (`predicted_rt`/`delta_rt_model`, populated later in
            // `retention_model::predict`) — see
            // `plans/lda_external_rt_iim_features.md`. 0.0 when
            // `--predicted-rt`/`--predicted-iim` aren't configured; the LDA
            // only includes the z² columns when they are (`ml/linear_discriminant.rs`).
            let (predicted_rt_external, delta_rt_z2_external, predicted_ims_external, delta_ims_z2_external) =
                self.external_z2(query, score.peptide.0 as usize, score.precursor_charge);

            // let (num_proteins, proteins) = self.db.assign_proteins(peptide);

            features.push(Feature {
                // Identifiers
                psm_id,
                peptide_idx: score.peptide,
                spec_id: query.id.clone(),
                file_id: query.file_id,
                rank: idx as u32 + 1,
                label: peptide.label(),
                expmass: precursor_mass,
                calcmass: peptide.monoisotopic,
                // Features
                charge: score.precursor_charge,
                rt: query.scan_start_time,
                ims: query
                    .precursors
                    .first()
                    .unwrap()
                    .inverse_ion_mobility
                    .unwrap_or(0.0),
                delta_mass,
                isotope_error,
                average_ppm: score.ppm_difference,
                hyperscore: score.hyperscore,
                delta_next: score.rank_key(self.ranking_score) - next,
                delta_best: best - score.rank_key(self.ranking_score),
                matched_peaks: k as u32,
                matched_intensity_pct: 100.0 * (score.summed_b + score.summed_y)
                    / query.total_ion_current,
                poisson: if log10_poisson.is_finite() {
                    log10_poisson
                } else {
                    f64::NEG_INFINITY
                },
                longest_b: score.longest_b as u32,
                longest_y: score.longest_y as u32,
                longest_y_pct: score.longest_y as f32 / (peptide.sequence.len() as f32),
                peptide_len: peptide.sequence.len(),
                scored_candidates: hits.scored_candidates as u32,
                missed_cleavages: peptide.missed_cleavages,

                // Outputs
                discriminant_score: 0.0,
                posterior_error: 1.0,
                spectrum_q: 1.0,
                protein_q: 1.0,
                peptide_q: 1.0,
                predicted_rt: 0.0,
                predicted_ims: 0.0,
                aligned_rt: query.scan_start_time,
                delta_rt_model: 0.999,
                delta_ims_model: 0.999,
                predicted_rt_external,
                delta_rt_z2_external,
                predicted_ims_external,
                delta_ims_z2_external,
                ms2_intensity: score.summed_b + score.summed_y,
                ms2_entropy_similarity: score.ms2_entropy_similarity,
                ms2_weighted_entropy_similarity: score.ms2_weighted_entropy_similarity,
                ms2_heuristic_entropy_similarity: score.ms2_heuristic_entropy_similarity,
                ms2_cosine_similarity: score.ms2_cosine_similarity,
                ms2_dot_product: score.ms2_dot_product,
                ms2_spectral_contrast_angle: score.ms2_spectral_contrast_angle,
                ms2_euclidean_similarity: score.ms2_euclidean_similarity,
                ms2_bray_curtis_similarity: score.ms2_bray_curtis_similarity,
                ms2_pearson_corr: score.ms2_pearson_corr,
                ms2_spearman_corr: score.ms2_spearman_corr,
                ms2_hypergeometric_probability: score.ms2_hypergeometric_probability,
                ms2_intersection: score.ms2_intersection,
                ms2_top6_matched_intensity: score.ms2_top6_matched_intensity,

                //Fragments
                protein_groups: None,
                num_protein_groups: 0,
                fragments,
                protein_group_q: 1.0,
            })
        }
    }

    /// Remove peaks matching a PSM from a query spectrum
    fn remove_matched_peaks(&self, query: &mut ProcessedSpectrum<Peak>, psm: &Feature) {
        let peptide = &self.db[psm.peptide_idx];
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind));

        let max_fragment_charge = max_fragment_charge(self.max_fragment_charge, psm.charge);

        // Remove MS2 peaks matched by previous match
        let mut to_remove = Vec::new();
        for frag in fragments {
            for charge in 1..max_fragment_charge {
                // Experimental peaks are multipled by charge, therefore theoretical are divided
                let mz = frag.monoisotopic_mass / charge as f32;
                if let Some(i) = crate::spectrum::select_most_intense_peak(
                    &query.peaks,
                    mz,
                    self.fragment_tol,
                    None,
                )
                {
                    to_remove.push(query.peaks[i]);
                }
            }
        }

        query.peaks = query
            .peaks
            .drain(..)
            .filter(|peak| !to_remove.contains(peak))
            .collect();
        query.total_ion_current = query.peaks.iter().map(|peak| peak.intensity).sum::<f32>();
    }

    /// Return multiple PSMs for each spectra - first is the best match, second PSM is the best match
    /// after all theoretical peaks assigned to the best match are removed, etc
    pub fn score_chimera_fast(&self, query: &ProcessedSpectrum<Peak>) -> Vec<Feature> {
        let precursor = query.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", query.id);
        });

        let mut query = query.clone();
        let hits = self.initial_hits(&query, precursor);

        let mut candidates: Vec<Feature> = Vec::with_capacity(self.report_psms);

        let mut prev = 0;
        while candidates.len() < self.report_psms {
            self.build_features(&query, precursor, &hits, 1, &mut candidates);
            if candidates.len() > prev {
                if let Some(feat) = candidates.get_mut(prev) {
                    self.remove_matched_peaks(&mut query, feat);
                    feat.rank = prev as u32 + 1;
                }
                prev = candidates.len()
            } else {
                break;
            }
        }
        candidates
    }

    /// Calculate full hyperscore for a given PSM
    fn score_candidate(
        &self,
        query: &ProcessedSpectrum<Peak>,
        pre_score: &PreScore,
    ) -> (Score, Option<Fragments>) {
        let mut score = Score {
            peptide: pre_score.peptide,
            precursor_charge: pre_score.precursor_charge,
            isotope_error: pre_score.isotope_error,
            ..Default::default()
        };
        let peptide = &self.db[score.peptide];
        let max_fragment_charge =
            max_fragment_charge(self.max_fragment_charge, score.precursor_charge);

        // Regenerate theoretical ions - initial database search might be
        // using only a subset of all possible ions (e.g. no b1/b2/y1/y2)
        // so we need to completely re-score this candidate
        let fragments = self
            .db
            .ion_kinds
            .iter()
            .flat_map(|kind| IonSeries::new(peptide, *kind).enumerate());

        let mut b_run = Run::default();
        let mut y_run = Run::default();

        let mut fragments_details = Fragments::default();

        // Unpack this candidate's sparse predicted MS2 vector (if
        // configured and this (peptide, charge) has a cache entry) into a
        // dense, zero-initialized `[f32; N_FRAGMENT_SLOTS]` once, up front
        // -- same "unpack sparse into a dense per-candidate buffer" shape
        // as the rest of this fork's dense-array precedent (see
        // `iim_dense_slot`). `observed_dense` is only allocated at all when
        // `predicted_dense` is `Some` -- no point tracking observed
        // intensities for a candidate with nothing to compare them against.
        // See `docs/ai/predicted_fragment_intensity.md`.
        let predicted_dense: Option<[f32; N_FRAGMENT_SLOTS]> = self
            .predicted_fragment_intensity_index
            .and_then(|by_slot| {
                let slot = Self::iim_dense_slot(
                    score.peptide.0 as usize,
                    score.precursor_charge,
                    self.min_precursor_charge,
                    self.max_precursor_charge,
                )?;
                by_slot.get(slot).copied().flatten()
            })
            .and_then(|(start, end)| {
                let annotation_id = self.predicted_fragment_intensity_annotation_id?;
                let intensity = self.predicted_fragment_intensity?;
                let (start, end) = (start as usize, end as usize);
                if start > end || end > annotation_id.len() || end > intensity.len() {
                    return None;
                }
                let mut dense = [0f32; N_FRAGMENT_SLOTS];
                for i in start..end {
                    if let Some(slot) = dense.get_mut(annotation_id[i] as usize) {
                        *slot = intensity[i].to_f32();
                    }
                }
                Some(dense)
            });
        let mut observed_dense = predicted_dense.map(|_| [0f32; N_FRAGMENT_SLOTS]);
        // Which of the 174 slots are a structurally real fragment position
        // for *this* peptide/charge -- set unconditionally per (idx, charge)
        // below, regardless of whether a peak matched there. Needed for the
        // sample-size-sensitive metrics (Pearson/Spearman correlation,
        // hypergeometric probability, intersection -- see
        // `docs/ai/predicted_fragment_intensity.md`), which must not treat
        // an unvisited slot the same as a real "no intensity" fragment.
        let mut is_real_dense = predicted_dense.map(|_| [false; N_FRAGMENT_SLOTS]);

        for (idx, frag) in fragments {
            for charge in 1..max_fragment_charge {
                let annotation_slot = ms2_similarity::fragment_annotation_id(frag.kind, idx, charge);
                if let (Some(is_real), Some(slot)) = (is_real_dense.as_mut(), annotation_slot) {
                    is_real[slot] = true;
                }

                // Experimental peaks are multipled by charge, therefore theoretical are divided
                let mz = frag.monoisotopic_mass / charge as f32;

                if let Some(i) = crate::spectrum::select_most_intense_peak(
                    &query.peaks,
                    mz,
                    self.fragment_tol,
                    None,
                ) {
                    let peak = &query.peaks[i];
                    let peak_charge = query.peak_charges.get(i).copied().unwrap_or(1);

                    score.ppm_difference +=
                        peak.intensity * (mz - peak.mass).abs() * 2E6 / (mz + peak.mass);

                    let exp_mz = peak.mass / peak_charge as f32 + PROTON;
                    let calc_mz = frag.monoisotopic_mass / peak_charge as f32 + PROTON;

                    match frag.kind {
                        Kind::A | Kind::B | Kind::C => {
                            score.matched_b += 1;
                            score.summed_b += peak.intensity;
                            b_run.matched(idx);
                        }
                        Kind::X | Kind::Y | Kind::Z => {
                            score.matched_y += 1;
                            score.summed_y += peak.intensity;
                            y_run.matched(idx);
                        }
                    }

                    if let (Some(observed_dense), Some(slot)) =
                        (observed_dense.as_mut(), annotation_slot)
                    {
                        observed_dense[slot] = peak.intensity;
                    }

                    if self.annotate_matches {
                        let idx = match frag.kind {
                            Kind::A | Kind::B | Kind::C => idx as i32 + 1,
                            Kind::X | Kind::Y | Kind::Z => {
                                peptide.sequence.len().saturating_sub(1) as i32 - idx as i32
                            }
                        };
                        fragments_details.kinds.push(frag.kind);
                        fragments_details.charges.push(peak_charge as i32);
                        fragments_details.mz_experimental.push(exp_mz);
                        fragments_details.mz_calculated.push(calc_mz);
                        fragments_details.fragment_ordinals.push(idx);
                        fragments_details.intensities.push(peak.intensity);
                    }
                }
            }
        }

        score.hyperscore = score.hyperscore(self.score_type);
        score.longest_b = b_run.longest;
        score.longest_y = y_run.longest;
        score.ppm_difference /= score.summed_b + score.summed_y;
        if let (Some(predicted), Some(observed), Some(is_real)) =
            (predicted_dense, observed_dense, is_real_dense)
        {
            // Compact once, feed every metric the same real-positions-only
            // pair -- not just Pearson/Spearman. The cache's `predicted`
            // covers all 3 Prosit fragment charges regardless of this
            // job's own `max_fragment_charge` (e.g. real F9477 production
            // config uses `max_fragment_charge: 1`, so SAGE's own loop
            // above only ever sets `is_real` for ~1/3 of the 174 slots) --
            // summing the *full* dense arrays would silently compare real
            // predicted intensities against phantom `observed = 0` for
            // fragment charges this job never even attempts to match,
            // biasing every metric downward. Found 2026-08-31, see
            // `docs/ai/predicted_fragment_intensity.md` and
            // `ms2_similarity::entropy_similarity`'s doc comment.
            let mut observed_real = Vec::with_capacity(N_FRAGMENT_SLOTS);
            let mut predicted_real = Vec::with_capacity(N_FRAGMENT_SLOTS);
            for i in 0..N_FRAGMENT_SLOTS {
                if is_real[i] {
                    observed_real.push(observed[i]);
                    predicted_real.push(predicted[i]);
                }
            }

            score.ms2_entropy_similarity =
                ms2_similarity::entropy_similarity(&observed_real, &predicted_real);
            score.ms2_weighted_entropy_similarity =
                ms2_similarity::weighted_entropy_similarity(&observed_real, &predicted_real);
            score.ms2_heuristic_entropy_similarity =
                ms2_similarity::heuristic_entropy_similarity(&observed_real, &predicted_real);
            score.ms2_cosine_similarity =
                ms2_similarity::cosine_similarity(&observed_real, &predicted_real);
            score.ms2_dot_product = ms2_similarity::dot_product(&observed_real, &predicted_real);
            score.ms2_spectral_contrast_angle =
                ms2_similarity::spectral_contrast_angle(&observed_real, &predicted_real);
            score.ms2_euclidean_similarity =
                ms2_similarity::euclidean_similarity(&observed_real, &predicted_real);
            score.ms2_bray_curtis_similarity =
                ms2_similarity::bray_curtis_similarity(&observed_real, &predicted_real);
            score.ms2_hypergeometric_probability =
                ms2_similarity::hypergeometric_probability(&observed_real, &predicted_real);
            score.ms2_intersection =
                ms2_similarity::intersection(&observed_real, &predicted_real, 20);
            let all_peak_intensities: Vec<f32> = query.peaks.iter().map(|p| p.intensity).collect();
            score.ms2_top6_matched_intensity = ms2_similarity::top6_matched_intensity(
                &observed_real,
                &predicted_real,
                &all_peak_intensities,
            );
            score.ms2_pearson_corr = ms2_similarity::pearson_corr(&observed_real, &predicted_real);
            score.ms2_spearman_corr = ms2_similarity::spearman_corr(&observed_real, &predicted_real);
        }

        if self.annotate_matches {
            (score, Some(fragments_details))
        } else {
            // drop(fragments_details);
            (score, None)
        }
    }
}

/// Maintain information about the longest continous ion ladder for a series
#[derive(Default)]
struct Run {
    start: usize,
    length: usize,
    last: usize,
    pub longest: usize,
}

impl Run {
    pub fn matched(&mut self, index: usize) {
        if self.last == index {
            return;
        } else if self.start + self.length == index {
            self.length += 1;
            self.longest = self.longest.max(self.length);
        } else {
            self.start = index;
            self.length = 1;
            self.longest = self.longest.max(self.length);
        }
        self.last = index;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn iim_dense_slot_in_range() {
        // charges 2..=4, charge_span 3 -- peptide 0's slots are 0,1,2 for
        // charge 2,3,4; peptide 1's are 3,4,5; etc.
        assert_eq!(Scorer::iim_dense_slot(0, 2, 2, 4), Some(0));
        assert_eq!(Scorer::iim_dense_slot(0, 3, 2, 4), Some(1));
        assert_eq!(Scorer::iim_dense_slot(0, 4, 2, 4), Some(2));
        assert_eq!(Scorer::iim_dense_slot(1, 2, 2, 4), Some(3));
        assert_eq!(Scorer::iim_dense_slot(5, 4, 2, 4), Some(17));
    }

    #[test]
    fn iim_dense_slot_out_of_range() {
        assert_eq!(Scorer::iim_dense_slot(0, 1, 2, 4), None);
        assert_eq!(Scorer::iim_dense_slot(0, 5, 2, 4), None);
    }

    #[test]
    fn iim_dense_slot_single_charge_span() {
        // charge_span 1 -- exactly the mk_scorer test fixture shape used
        // in crates/sage/tests/integration.rs's IIM eviction tests.
        assert_eq!(Scorer::iim_dense_slot(0, 1, 1, 1), Some(0));
        assert_eq!(Scorer::iim_dense_slot(7, 1, 1, 1), Some(7));
        assert_eq!(Scorer::iim_dense_slot(0, 2, 1, 1), None);
    }

    #[test]
    fn longest_series() {
        let mut run = Run::default();

        run.matched(1);
        run.matched(2);
        run.matched(3);
        run.matched(3);
        run.matched(3);

        assert_eq!(run.length, 3);
        assert_eq!(run.longest, 3);

        run.matched(5);
        run.matched(5);
        assert_eq!(run.length, 1);
        assert_eq!(run.longest, 3);
        run.matched(6);
        assert_eq!(run.length, 2);
    }

    #[test]
    fn test_max_fragment_charge() {
        assert_eq!(max_fragment_charge(None, 1), 2);
        assert_eq!(max_fragment_charge(None, 2), 2);
        assert_eq!(max_fragment_charge(None, 3), 3);
        assert_eq!(max_fragment_charge(None, 4), 4);
        assert_eq!(max_fragment_charge(Some(1), 2), 2);
        assert_eq!(max_fragment_charge(Some(1), 3), 2);
        assert_eq!(max_fragment_charge(Some(2), 4), 3);
        assert_eq!(max_fragment_charge(Some(4), 1), 2);
    }
}
