//! MS2 fragment-intensity similarity, computed from `git/featureprediction`'s
//! predicted intensities against observed matched-peak intensities -- see
//! `docs/ai/predicted_fragment_intensity.md`. Feature-only: nothing here
//! evicts or re-ranks candidates (contrast `evict_rt_iim_mismatches`/
//! `Score::combined_score` in `scoring.rs`).

use crate::ion_series::Kind;

/// Size of the dense per-candidate predicted/observed intensity vectors --
/// matches `git/featureprediction`'s `fragment_slots.N_SLOTS` exactly: 29
/// ordinals x 2 kinds (b, y) x 3 fragment charges.
pub const N_FRAGMENT_SLOTS: usize = 174;

const N_ORDINALS: usize = 29;
const N_CHARGES: u8 = 3;

/// `kind_block(kind) + idx*3 + (charge-1)` -- this repo's SAGE-order
/// fragment-slot numbering (see `git/featureprediction`'s
/// `fragment_slots.py::annotation_id_for`, which this must stay bit-for-bit
/// consistent with). `idx` is SAGE's own raw backbone-position counter from
/// `IonSeries::enumerate()` (0-based, shared between b/y series) -- exactly
/// what `score_candidate`'s existing fragment loop already has, no ordinal
/// recomputation needed. `None` for any fragment outside the Prosit
/// `compact_trt` model's vocabulary (kinds other than b/y, `idx >= 29`,
/// fragment charge `> 3`) -- those fragments simply have no predicted
/// intensity, same treatment as a missing cache entry.
pub fn fragment_annotation_id(kind: Kind, idx: usize, charge: u8) -> Option<usize> {
    if idx >= N_ORDINALS || charge == 0 || charge > N_CHARGES {
        return None;
    }
    let kind_block = match kind {
        Kind::B => 0,
        Kind::Y => N_ORDINALS * N_CHARGES as usize,
        _ => return None,
    };
    Some(kind_block + idx * N_CHARGES as usize + (charge - 1) as usize)
}

fn shannon_entropy_term(p: f32) -> f32 {
    if p > 0.0 {
        -p * p.ln()
    } else {
        0.0
    }
}

/// Unweighted spectral entropy similarity (Li et al. 2021, "Spectral
/// entropy outperforms MS/MS dot product similarity for small-molecule
/// compound identification", Nature Methods) between two intensity vectors
/// addressed by the same fixed vocabulary (this repo's `fragment_annotation_id`
/// slots) -- `1 - JSD(P, Q) / ln(2)`, where `P`/`Q` are `observed`/`predicted`
/// each normalized to sum to 1, and `JSD` is their Jensen-Shannon divergence
/// (`H((P+Q)/2) - 0.5*H(P) - 0.5*H(Q)`). Ranges `[0, 1]`: `1.0` for
/// identical (post-normalization) distributions, `0.0` for disjoint support
/// or either vector being all-zero (empty spectra have no meaningful
/// similarity, not perfect similarity).
///
/// **Both slices must already be restricted to this candidate's real
/// fragment positions** (`score_candidate`'s `is_real_dense` mask, applied
/// once via a single compaction pass shared by every metric in this
/// module) -- **not** the full `N_FRAGMENT_SLOTS` dense array. This isn't
/// just an efficiency choice: the cache's `predicted` values cover all 3
/// Prosit fragment charges regardless of what this job's own
/// `max_fragment_charge` config actually searches (e.g. the real F9477
/// production config uses `max_fragment_charge: 1`, so SAGE's own fragment
/// loop only ever visits ~1/3 of the 174-slot vocabulary) -- summing over
/// the full unfiltered array would silently compare real predicted
/// intensities against phantom `observed = 0` values for fragment charges
/// this job never even attempts to match, biasing every metric downward.
/// Found and fixed 2026-08-31 (see `docs/ai/predicted_fragment_intensity.md`).
pub fn entropy_similarity(observed: &[f32], predicted: &[f32]) -> f32 {
    let sum_observed: f32 = observed.iter().sum();
    let sum_predicted: f32 = predicted.iter().sum();
    if sum_observed <= 0.0 || sum_predicted <= 0.0 {
        return 0.0;
    }

    let mut h_p = 0.0f32;
    let mut h_q = 0.0f32;
    let mut h_m = 0.0f32;
    for i in 0..observed.len() {
        let p = observed[i] / sum_observed;
        let q = predicted[i] / sum_predicted;
        let m = 0.5 * (p + q);
        h_p += shannon_entropy_term(p);
        h_q += shannon_entropy_term(q);
        h_m += shannon_entropy_term(m);
    }
    let jsd = h_m - 0.5 * h_p - 0.5 * h_q;
    (1.0 - jsd / std::f32::consts::LN_2).clamp(0.0, 1.0)
}

fn unit_normalize(v: &[f32]) -> Vec<f32> {
    let sum_sq: f32 = v.iter().map(|x| x * x).sum();
    if sum_sq <= 0.0 {
        return vec![0.0f32; v.len()];
    }
    let norm = sum_sq.sqrt();
    v.iter().map(|x| x / norm).collect()
}

/// `H(one_normalize(v))` without materializing the normalized copy --
/// `-Σ (v[i]/sum)·ln(v[i]/sum)`.
fn shannon_entropy_one_normalized(v: &[f32]) -> f32 {
    let sum: f32 = v.iter().sum();
    if sum <= 0.0 {
        return 0.0;
    }
    v.iter().map(|&x| shannon_entropy_term(x / sum)).sum()
}

/// MSBooster `SpectrumComparison::cosineSimilarity` -- raw (non-normalized)
/// vectors: `Σ(observed·predicted) / sqrt(Σobserved² · Σpredicted²)`. `0.0`
/// if either sum-of-squares is zero (undefined direction). See
/// [`entropy_similarity`]'s doc comment for why both slices must already be
/// restricted to real fragment positions.
pub fn cosine_similarity(observed: &[f32], predicted: &[f32]) -> f32 {
    let mut dot = 0.0f32;
    let mut sum_observed_sq = 0.0f32;
    let mut sum_predicted_sq = 0.0f32;
    for i in 0..observed.len() {
        dot += observed[i] * predicted[i];
        sum_observed_sq += observed[i] * observed[i];
        sum_predicted_sq += predicted[i] * predicted[i];
    }
    let denom = (sum_observed_sq * sum_predicted_sq).sqrt();
    if denom > 0.0 {
        dot / denom
    } else {
        0.0
    }
}

/// MSBooster `SpectrumComparison::dotProduct` -- algebraically the cosine
/// similarity of the two L2-normalized vectors (so numerically identical to
/// [`cosine_similarity`] whenever both are well-defined); kept as a
/// separate function to preserve MSBooster's own zero-check nuance: `0.0`
/// if `observed`'s raw (pre-normalization) sum-of-squares is zero, checked
/// independently of [`cosine_similarity`]'s own (symmetric) zero-check.
pub fn dot_product(observed: &[f32], predicted: &[f32]) -> f32 {
    let sum_observed_sq: f32 = observed.iter().map(|x| x * x).sum();
    if sum_observed_sq <= 0.0 {
        return 0.0;
    }
    let unit_observed = unit_normalize(observed);
    let unit_predicted = unit_normalize(predicted);
    unit_observed
        .iter()
        .zip(unit_predicted.iter())
        .map(|(a, b)| a * b)
        .sum()
}

/// MSBooster `SpectrumComparison::spectralContrastAngle` -- angular form of
/// [`cosine_similarity`]: `1 - (2/π)·acos(cosine)`. Cosine is clamped to
/// `[-1, 1]` before `acos` -- floating-point summation error can otherwise
/// push it fractionally outside that range and produce `NaN`.
pub fn spectral_contrast_angle(observed: &[f32], predicted: &[f32]) -> f32 {
    let cosine = cosine_similarity(observed, predicted).clamp(-1.0, 1.0);
    1.0 - (2.0 / std::f32::consts::PI) * cosine.acos()
}

/// MSBooster `SpectrumComparison::euclideanDistance` -- despite the name,
/// returns a similarity-shaped value: `1 - ||unit(observed) -
/// unit(predicted)||₂` on L2-normalized vectors. `1 - sqrt(2)` (the
/// maximum possible distance between two unit vectors) if `observed`'s raw
/// sum-of-squares is zero.
pub fn euclidean_similarity(observed: &[f32], predicted: &[f32]) -> f32 {
    let sum_observed_sq: f32 = observed.iter().map(|x| x * x).sum();
    if sum_observed_sq <= 0.0 {
        return 1.0 - std::f32::consts::SQRT_2;
    }
    let unit_observed = unit_normalize(observed);
    let unit_predicted = unit_normalize(predicted);
    let sq_dist: f32 = unit_observed
        .iter()
        .zip(unit_predicted.iter())
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    1.0 - sq_dist.sqrt()
}

/// MSBooster `SpectrumComparison::brayCurtis` -- similarity-shaped
/// (`1 - dissimilarity`) on L2-normalized vectors: `1 - (Σ|o-p| /
/// Σ(o+p))`. `0.0` (not `1.0`) if `observed`'s raw sum-of-squares is zero.
pub fn bray_curtis_similarity(observed: &[f32], predicted: &[f32]) -> f32 {
    let sum_observed_sq: f32 = observed.iter().map(|x| x * x).sum();
    if sum_observed_sq <= 0.0 {
        return 0.0;
    }
    let unit_observed = unit_normalize(observed);
    let unit_predicted = unit_normalize(predicted);
    let mut num = 0.0f32;
    let mut den = 0.0f32;
    for i in 0..unit_observed.len() {
        num += (unit_observed[i] - unit_predicted[i]).abs();
        den += unit_observed[i] + unit_predicted[i];
    }
    if den > 0.0 {
        1.0 - num / den
    } else {
        0.0
    }
}

/// MSBooster `SpectrumComparison::weightedSpectralEntropy` -- **not**
/// weighted by the m/z-frequency table the `weighted*` metrics above use
/// (that weighting source, `getWeights`/`freqs`, is unverified -- see
/// `docs/ai/predicted_fragment_intensity.md` -- so those metrics are
/// skipped entirely for now). This is a self-weighting transform using the
/// predicted spectrum's own entropy: [`entropy_similarity`] multiplied by
/// `H(one_normalize(predicted))^0.5`. Not bounded to `[0, 1]` -- the
/// multiplier can exceed `1.0` for a high-entropy (flat/noisy) predicted
/// spectrum, matching MSBooster's own unbounded behavior.
pub fn weighted_entropy_similarity(observed: &[f32], predicted: &[f32]) -> f32 {
    let base = entropy_similarity(observed, predicted);
    if base <= 0.0 {
        return 0.0;
    }
    let h_predicted = shannon_entropy_one_normalized(predicted);
    base * h_predicted.sqrt()
}

/// MSBooster `SpectrumComparison::heuristicSpectralEntropy` -- if the
/// predicted spectrum's own entropy is low (`< 1.75`, i.e. concentrated on
/// few fragments), both vectors are reweighted element-wise by
/// `intensity^(H/2.75)` before recomputing [`entropy_similarity`];
/// otherwise falls back to plain [`entropy_similarity`] unchanged. `0.0`
/// intensities are left at `0.0` rather than `powf`'d (`0f32.powf(0.0) ==
/// 1.0` in IEEE754, which would wrongly turn a true absence into a fake
/// positive value at the `H == 0` boundary).
pub fn heuristic_entropy_similarity(observed: &[f32], predicted: &[f32]) -> f32 {
    let h_predicted = shannon_entropy_one_normalized(predicted);
    if h_predicted >= 1.75 {
        return entropy_similarity(observed, predicted);
    }
    let exponent = h_predicted / 2.75;
    let reweight = |x: f32| if x > 0.0 { x.powf(exponent) } else { 0.0 };
    let reweighted_observed: Vec<f32> = observed.iter().copied().map(reweight).collect();
    let reweighted_predicted: Vec<f32> = predicted.iter().copied().map(reweight).collect();
    entropy_similarity(&reweighted_observed, &reweighted_predicted)
}

/// MSBooster `SpectrumComparison::pearsonCorr` -- plain Pearson correlation
/// over two equal-length slices (the caller restricts these to this
/// candidate's *real* fragment positions only, not the full
/// `N_FRAGMENT_SLOTS` dense array -- Pearson/Spearman are sample-size
/// sensitive, unlike every metric above, so padding with extra `(0, 0)`
/// slots would silently bias the result, not just add neutral terms). `-1.0`
/// (minimum-correlation sentinel, matching MSBooster's own convention) if
/// fewer than 2 positions or either side has zero variance.
pub fn pearson_corr(observed: &[f32], predicted: &[f32]) -> f32 {
    let n = observed.len();
    if n < 2 || n != predicted.len() {
        return -1.0;
    }
    let mean_observed = observed.iter().sum::<f32>() / n as f32;
    let mean_predicted = predicted.iter().sum::<f32>() / n as f32;
    let mut cov = 0.0f32;
    let mut var_observed = 0.0f32;
    let mut var_predicted = 0.0f32;
    for i in 0..n {
        let d_observed = observed[i] - mean_observed;
        let d_predicted = predicted[i] - mean_predicted;
        cov += d_observed * d_predicted;
        var_observed += d_observed * d_observed;
        var_predicted += d_predicted * d_predicted;
    }
    let denom = (var_observed * var_predicted).sqrt();
    if denom > 0.0 {
        cov / denom
    } else {
        -1.0
    }
}

/// 1-based ranks with average-rank tie handling (matching Apache Commons
/// `NaturalRanking`'s default strategy, which MSBooster's `spearmanCorr`
/// uses via `SpearmansCorrelation`).
fn rank_with_average_ties(v: &[f32]) -> Vec<f32> {
    let n = v.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut ranks = vec![0.0f32; n];
    let mut i = 0;
    while i < n {
        let mut j = i;
        while j + 1 < n && v[order[j + 1]] == v[order[i]] {
            j += 1;
        }
        let average_rank = ((i + 1) + (j + 1)) as f32 / 2.0;
        for slot in order.iter().take(j + 1).skip(i) {
            ranks[*slot] = average_rank;
        }
        i = j + 1;
    }
    ranks
}

/// MSBooster `SpectrumComparison::spearmanCorr` -- [`pearson_corr`] of the
/// two vectors' ranks. Same real-positions-only restriction and `-1.0`
/// sentinel as `pearson_corr`.
pub fn spearman_corr(observed: &[f32], predicted: &[f32]) -> f32 {
    if observed.len() < 2 || observed.len() != predicted.len() {
        return -1.0;
    }
    pearson_corr(
        &rank_with_average_ties(observed),
        &rank_with_average_ties(predicted),
    )
}

fn ln_factorial(n: u64) -> f64 {
    (2..=n).map(|i| (i as f64).ln()).sum()
}

fn ln_choose(n: u64, k: u64) -> f64 {
    if k > n {
        f64::NEG_INFINITY
    } else {
        ln_factorial(n) - ln_factorial(k) - ln_factorial(n - k)
    }
}

/// `P(X >= successes)` for `X ~ Hypergeometric(population, population_successes,
/// sample_size)`, via direct log-space summation of the PMF -- population
/// sizes here are always small (bounded by this repo's 174-slot vocabulary,
/// see `hypergeometric_probability` below), so an exact `ln_factorial` loop
/// is both cheap and more accurate than a Stirling approximation at this
/// scale, and no external stats crate is needed.
fn hypergeometric_upper_tail(
    population: u64,
    population_successes: u64,
    sample_size: u64,
    successes: u64,
) -> f64 {
    if population == 0 || sample_size == 0 || sample_size > population {
        return 1.0;
    }
    let population_failures = population - population_successes;
    let lo = successes.max(sample_size.saturating_sub(population_failures));
    let hi = sample_size.min(population_successes);
    if lo > hi {
        return 0.0;
    }
    let ln_denom = ln_choose(population, sample_size);
    let p: f64 = (lo..=hi)
        .map(|i| (ln_choose(population_successes, i) + ln_choose(population_failures, sample_size - i) - ln_denom).exp())
        .sum();
    p.min(1.0)
}

/// MSBooster `SpectrumComparison::hypergeometricProbability`, **restricted
/// in scope** (per explicit decision, see
/// `docs/ai/predicted_fragment_intensity.md`): MSBooster's "population" is
/// every theoretically possible fragment across all 6 ion kinds within the
/// scan's observable m/z range; this instead uses this repo's own 174-slot
/// vocabulary, restricted to this candidate's real positions (`observed`/
/// `predicted` must already be compacted to just those, same as every
/// other metric in this module -- see [`entropy_similarity`]'s doc
/// comment) -- narrower than MSBooster's literal definition, but avoids
/// generating fragment kinds beyond what the job's own search already
/// uses. Still a coherent, non-degenerate test: "population" = real
/// positions (i.e. `observed.len()`), "population successes" = how many
/// were observed at all, "sample" = the subset the prediction model
/// actually covered (nonzero `predicted`), "sample successes" = how many
/// of *those* were observed. Returns `-log10(P(X >= sample_successes))`,
/// `0.0` if population or sample is empty.
pub fn hypergeometric_probability(observed: &[f32], predicted: &[f32]) -> f32 {
    let population = observed.len() as u64;
    if population == 0 {
        return 0.0;
    }
    let mut population_successes = 0u64;
    let mut sample_size = 0u64;
    let mut sample_successes = 0u64;
    for i in 0..observed.len() {
        let matched = observed[i] > 0.0;
        if matched {
            population_successes += 1;
        }
        if predicted[i] > 0.0 {
            sample_size += 1;
            if matched {
                sample_successes += 1;
            }
        }
    }
    if sample_size == 0 {
        return 0.0;
    }
    let p = hypergeometric_upper_tail(population, population_successes, sample_size, sample_successes)
        .max(1e-300);
    (-p.log10()) as f32
}

/// MSBooster `SpectrumComparison::intersection`, **restricted in scope**
/// (same reasoning as [`hypergeometric_probability`]): population = this
/// candidate's real positions within this repo's 174-slot vocabulary (see
/// [`entropy_similarity`]'s doc comment -- `observed`/`predicted` must
/// already be compacted to those), not MSBooster's literal "all 6 ion
/// kinds" set. Sorts `observed` descending, takes the top `top_n`, counts
/// how many also have `predicted > 0`. Index-aligned membership, not
/// MSBooster's raw float-value-membership check -- this repo's arrays are
/// already index-aligned by construction, so value-equality matching
/// (needed in MSBooster's own less-aligned data structures) isn't
/// necessary here. Returns a raw count, not a fraction (matching
/// MSBooster: explicitly not Jaccard).
pub fn intersection(observed: &[f32], predicted: &[f32], top_n: usize) -> u32 {
    let mut indexed: Vec<(usize, f32)> = observed.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
    indexed
        .iter()
        .take(top_n)
        .filter(|&&(i, _)| predicted[i] > 0.0)
        .count() as u32
}

/// MSBooster `SpectrumComparison::top6matchedIntensity`, **restricted in
/// scope** for its numerator only (same reasoning as
/// [`hypergeometric_probability`]/[`intersection`] -- `observed`/
/// `predicted` must already be compacted to real positions): picks the
/// top 6 (or fewer) real positions by *predicted* intensity, sums
/// `ln(observed + 1)` over those. The denominator is **not**
/// scope-restricted -- it needs the full raw spectrum (`Σ ln(peak_intensity
/// + 1)` over every observed peak above the scan's own mean intensity,
/// matched-to-a-fragment or not), which `all_observed_peak_intensities`
/// (the caller's `query.peaks` intensities, unfiltered) provides directly.
/// `0.0` if there are no peaks at all, or the denominator isn't `> 0`
/// (covers a degenerate all-zero/uniform spectrum too).
pub fn top6_matched_intensity(
    observed: &[f32],
    predicted: &[f32],
    all_observed_peak_intensities: &[f32],
) -> f32 {
    let mut indexed: Vec<(usize, f32)> = predicted.iter().copied().enumerate().collect();
    indexed.sort_by(|a, b| b.1.total_cmp(&a.1));
    let numerator: f32 = indexed
        .iter()
        .take(6)
        .map(|&(i, _)| (observed[i] + 1.0).ln())
        .sum();

    if all_observed_peak_intensities.is_empty() {
        return 0.0;
    }
    let mean = all_observed_peak_intensities.iter().sum::<f32>()
        / all_observed_peak_intensities.len() as f32;
    let denominator: f32 = all_observed_peak_intensities
        .iter()
        .filter(|&&i| i > mean)
        .map(|&i| (i + 1.0).ln())
        .sum();
    if denominator > 0.0 {
        numerator / denominator
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn annotation_id_matches_python_kind_blocks() {
        // b_idx0+1 -- see fragment_slots.py's annotation_names()[0]
        assert_eq!(fragment_annotation_id(Kind::B, 0, 1), Some(0));
        // b_idx0+3
        assert_eq!(fragment_annotation_id(Kind::B, 0, 3), Some(2));
        // b_idx1+1
        assert_eq!(fragment_annotation_id(Kind::B, 1, 1), Some(3));
        // y block starts at 87
        assert_eq!(fragment_annotation_id(Kind::Y, 0, 1), Some(87));
        // last valid slot: y, idx=28, charge=3 -> 87 + 28*3 + 2 = 173
        assert_eq!(fragment_annotation_id(Kind::Y, 28, 3), Some(173));
    }

    #[test]
    fn annotation_id_out_of_vocabulary_is_none() {
        assert_eq!(fragment_annotation_id(Kind::B, 29, 1), None);
        assert_eq!(fragment_annotation_id(Kind::B, 0, 4), None);
        assert_eq!(fragment_annotation_id(Kind::B, 0, 0), None);
        assert_eq!(fragment_annotation_id(Kind::A, 0, 1), None);
        assert_eq!(fragment_annotation_id(Kind::X, 0, 1), None);
        assert_eq!(fragment_annotation_id(Kind::Z, 0, 1), None);
    }

    #[test]
    fn entropy_similarity_identical_vectors_is_one() {
        let mut v = [0.0f32; N_FRAGMENT_SLOTS];
        v[0] = 0.3;
        v[10] = 0.7;
        v[100] = 1.0;
        assert!((entropy_similarity(&v, &v) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn entropy_similarity_disjoint_support_is_zero() {
        let mut observed = [0.0f32; N_FRAGMENT_SLOTS];
        let mut predicted = [0.0f32; N_FRAGMENT_SLOTS];
        observed[0] = 1.0;
        predicted[1] = 1.0;
        assert!(entropy_similarity(&observed, &predicted).abs() < 1e-5);
    }

    #[test]
    fn entropy_similarity_is_scale_invariant() {
        let mut observed = [0.0f32; N_FRAGMENT_SLOTS];
        observed[0] = 1.0;
        observed[5] = 3.0;
        let mut predicted = observed;
        for v in predicted.iter_mut() {
            *v *= 1000.0;
        }
        let a = entropy_similarity(&observed, &predicted);
        assert!((a - 1.0).abs() < 1e-4);
    }

    #[test]
    fn entropy_similarity_all_zero_vector_is_zero_not_one() {
        let zeros = [0.0f32; N_FRAGMENT_SLOTS];
        let mut nonzero = [0.0f32; N_FRAGMENT_SLOTS];
        nonzero[0] = 1.0;
        assert_eq!(entropy_similarity(&zeros, &nonzero), 0.0);
        assert_eq!(entropy_similarity(&zeros, &zeros), 0.0);
    }

    #[test]
    fn entropy_similarity_partial_overlap_between_zero_and_one() {
        let mut observed = [0.0f32; N_FRAGMENT_SLOTS];
        let mut predicted = [0.0f32; N_FRAGMENT_SLOTS];
        observed[0] = 1.0;
        observed[1] = 1.0;
        predicted[0] = 1.0;
        predicted[2] = 1.0;
        let s = entropy_similarity(&observed, &predicted);
        assert!(s > 0.0 && s < 1.0, "expected partial similarity, got {s}");
    }

    fn make(entries: &[(usize, f32)]) -> [f32; N_FRAGMENT_SLOTS] {
        let mut v = [0.0f32; N_FRAGMENT_SLOTS];
        for &(i, val) in entries {
            v[i] = val;
        }
        v
    }

    #[test]
    fn cosine_similarity_identical_is_one() {
        let v = make(&[(0, 0.3), (10, 0.7), (100, 1.0)]);
        assert!((cosine_similarity(&v, &v) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn cosine_similarity_orthogonal_is_zero() {
        let a = make(&[(0, 1.0)]);
        let b = make(&[(1, 1.0)]);
        assert_eq!(cosine_similarity(&a, &b), 0.0);
    }

    #[test]
    fn cosine_similarity_scale_invariant() {
        let a = make(&[(0, 1.0), (5, 3.0)]);
        let mut b = a;
        for v in b.iter_mut() {
            *v *= 42.0;
        }
        assert!((cosine_similarity(&a, &b) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn cosine_similarity_zero_vector_is_zero() {
        let zeros = [0.0f32; N_FRAGMENT_SLOTS];
        let v = make(&[(0, 1.0)]);
        assert_eq!(cosine_similarity(&zeros, &v), 0.0);
    }

    #[test]
    fn dot_product_matches_cosine_when_well_defined() {
        let a = make(&[(0, 1.0), (5, 3.0)]);
        let b = make(&[(0, 2.0), (5, 1.0), (9, 4.0)]);
        assert!((dot_product(&a, &b) - cosine_similarity(&a, &b)).abs() < 1e-5);
    }

    #[test]
    fn dot_product_zero_observed_is_zero() {
        let zeros = [0.0f32; N_FRAGMENT_SLOTS];
        let v = make(&[(0, 1.0)]);
        assert_eq!(dot_product(&zeros, &v), 0.0);
    }

    #[test]
    fn spectral_contrast_angle_identical_is_one() {
        let v = make(&[(0, 1.0), (5, 2.0)]);
        assert!((spectral_contrast_angle(&v, &v) - 1.0).abs() < 1e-4);
    }

    #[test]
    fn spectral_contrast_angle_orthogonal_is_zero() {
        let a = make(&[(0, 1.0)]);
        let b = make(&[(1, 1.0)]);
        assert!(spectral_contrast_angle(&a, &b).abs() < 1e-4);
    }

    #[test]
    fn euclidean_similarity_identical_is_one() {
        let v = make(&[(0, 1.0), (5, 2.0)]);
        assert!((euclidean_similarity(&v, &v) - 1.0).abs() < 1e-4);
    }

    #[test]
    fn euclidean_similarity_zero_observed_is_min() {
        let zeros = [0.0f32; N_FRAGMENT_SLOTS];
        let v = make(&[(0, 1.0)]);
        let expected = 1.0 - std::f32::consts::SQRT_2;
        assert!((euclidean_similarity(&zeros, &v) - expected).abs() < 1e-5);
    }

    #[test]
    fn bray_curtis_similarity_identical_is_one() {
        let v = make(&[(0, 1.0), (5, 2.0)]);
        assert!((bray_curtis_similarity(&v, &v) - 1.0).abs() < 1e-4);
    }

    #[test]
    fn bray_curtis_similarity_disjoint_is_zero() {
        let a = make(&[(0, 1.0)]);
        let b = make(&[(1, 1.0)]);
        assert!(bray_curtis_similarity(&a, &b).abs() < 1e-4);
    }

    #[test]
    fn weighted_entropy_similarity_zero_when_base_is_zero() {
        let a = make(&[(0, 1.0)]);
        let b = make(&[(1, 1.0)]);
        assert_eq!(weighted_entropy_similarity(&a, &b), 0.0);
    }

    #[test]
    fn weighted_entropy_similarity_positive_for_partial_overlap() {
        let observed = make(&[(0, 1.0), (1, 1.0)]);
        let predicted = make(&[(0, 1.0), (2, 1.0)]);
        assert!(weighted_entropy_similarity(&observed, &predicted) > 0.0);
    }

    #[test]
    fn heuristic_entropy_similarity_matches_base_for_high_entropy_prediction() {
        // Flat, high-entropy predicted spectrum (>=1.75 nats) -- falls back
        // to plain entropy_similarity unchanged.
        let mut predicted = [0.0f32; N_FRAGMENT_SLOTS];
        for i in 0..20 {
            predicted[i] = 1.0;
        }
        let observed = predicted;
        assert_eq!(
            heuristic_entropy_similarity(&observed, &predicted),
            entropy_similarity(&observed, &predicted)
        );
    }

    #[test]
    fn heuristic_entropy_similarity_reweights_low_entropy_prediction() {
        // Concentrated (low-entropy) predicted spectrum -- reweight path
        // taken; identical vectors should still be perfectly self-similar.
        let v = make(&[(0, 1.0)]);
        assert!((heuristic_entropy_similarity(&v, &v) - 1.0).abs() < 1e-4);
    }

    #[test]
    fn pearson_corr_identical_is_one() {
        let v = [1.0f32, 2.0, 3.0, 4.0];
        assert!((pearson_corr(&v, &v) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn pearson_corr_inverted_is_negative_one() {
        let a = [1.0f32, 2.0, 3.0, 4.0];
        let b = [4.0f32, 3.0, 2.0, 1.0];
        assert!((pearson_corr(&a, &b) - (-1.0)).abs() < 1e-5);
    }

    #[test]
    fn pearson_corr_too_short_is_sentinel() {
        assert_eq!(pearson_corr(&[1.0], &[1.0]), -1.0);
        assert_eq!(pearson_corr(&[], &[]), -1.0);
    }

    #[test]
    fn pearson_corr_zero_variance_is_sentinel() {
        let a = [1.0f32, 1.0, 1.0];
        let b = [1.0f32, 2.0, 3.0];
        assert_eq!(pearson_corr(&a, &b), -1.0);
    }

    #[test]
    fn spearman_corr_identical_is_one() {
        let v = [1.0f32, 5.0, 2.0, 8.0];
        assert!((spearman_corr(&v, &v) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn spearman_corr_monotonic_nonlinear_is_one() {
        // Spearman only cares about rank order, not linear scale.
        let a = [1.0f32, 2.0, 3.0, 4.0];
        let b = [1.0f32, 10.0, 1000.0, 100000.0];
        assert!((spearman_corr(&a, &b) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn spearman_corr_handles_ties() {
        let a = [1.0f32, 1.0, 2.0, 3.0];
        let b = [1.0f32, 1.0, 2.0, 3.0];
        assert!((spearman_corr(&a, &b) - 1.0).abs() < 1e-5);
    }

    #[test]
    fn hypergeometric_probability_zero_without_population() {
        assert_eq!(hypergeometric_probability(&[], &[]), 0.0);
    }

    #[test]
    fn hypergeometric_probability_higher_for_more_surprising_match() {
        // 20 real positions; predicted covers 5 of them (the sample).
        // Case A: all 5 sample positions are observed-matched, plus 2 more
        // matches outside the sample (population has 7 total matches) --
        // a small sample matching almost perfectly is surprising.
        // Case B: same 7 total population matches, but only 1 of them
        // happens to fall inside the 5-slot sample -- much less surprising.
        let mut predicted = [0.0f32; 20];
        predicted[0..5].fill(1.0);

        let mut observed_surprising = [0.0f32; 20];
        observed_surprising[0..7].fill(1.0); // slots 0..5 (in sample) + 5,6 (outside)
        let mut observed_unsurprising = [0.0f32; 20];
        observed_unsurprising[0] = 1.0; // 1 of 5 sample slots
        observed_unsurprising[10..16].fill(1.0); // 6 more matches outside the sample

        let p_surprising = hypergeometric_probability(&observed_surprising, &predicted);
        let p_unsurprising = hypergeometric_probability(&observed_unsurprising, &predicted);
        assert!(
            p_surprising > p_unsurprising,
            "expected {p_surprising} > {p_unsurprising}"
        );
    }

    #[test]
    fn intersection_counts_top_overlap() {
        let observed = [5.0f32, 4.0, 3.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0];
        let mut predicted = [0.0f32; 10];
        predicted[0] = 1.0;
        predicted[2] = 1.0;
        predicted[9] = 1.0; // slot 9 not in top-3 observed
        assert_eq!(intersection(&observed, &predicted, 3), 2);
    }

    #[test]
    fn top6_matched_intensity_zero_without_peaks() {
        assert_eq!(top6_matched_intensity(&[], &[], &[]), 0.0);
    }

    #[test]
    fn top6_matched_intensity_higher_when_top_predicted_slots_are_observed() {
        // predicted's top-6 by value: slots 0..6 (value 6.0 down to 1.0).
        let predicted = [6.0f32, 5.0, 4.0, 3.0, 2.0, 1.0, 0.5, 0.5, 0.5, 0.5];
        let observed_matches_top6 = [9.0f32, 9.0, 9.0, 9.0, 9.0, 9.0, 0.0, 0.0, 0.0, 0.0];
        let observed_misses_top6 = [0.0f32, 0.0, 0.0, 0.0, 0.0, 0.0, 9.0, 9.0, 9.0, 9.0];
        let peaks = [9.0f32, 1.0, 1.0, 1.0]; // mean=3.0, only the 9.0 peak is "above mean"

        let higher = top6_matched_intensity(&observed_matches_top6, &predicted, &peaks);
        let lower = top6_matched_intensity(&observed_misses_top6, &predicted, &peaks);
        assert!(higher > lower, "expected {higher} > {lower}");
    }

    #[test]
    fn intersection_zero_without_overlap() {
        // top_n=1 restricts to only the single most-intense observed
        // position (slot 0); predicted has no value there, so no overlap
        // -- unlike top_n>=population, which would include every real
        // position regardless of its own observed value.
        let observed = [1.0f32, 0.0];
        let predicted = [0.0f32, 1.0];
        assert_eq!(intersection(&observed, &predicted, 1), 0);
    }
}
