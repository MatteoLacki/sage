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
/// A slot where both `observed` and `predicted` are `0.0` (e.g. a fragment
/// position beyond this peptide's real length, or outside the model's
/// vocabulary) contributes nothing to any of the three entropies -- safe to
/// leave both dense arrays zero-initialized for slots this peptide never
/// visits, no separate "valid slot" mask needed.
pub fn entropy_similarity(
    observed: &[f32; N_FRAGMENT_SLOTS],
    predicted: &[f32; N_FRAGMENT_SLOTS],
) -> f32 {
    let sum_observed: f32 = observed.iter().sum();
    let sum_predicted: f32 = predicted.iter().sum();
    if sum_observed <= 0.0 || sum_predicted <= 0.0 {
        return 0.0;
    }

    let mut h_p = 0.0f32;
    let mut h_q = 0.0f32;
    let mut h_m = 0.0f32;
    for i in 0..N_FRAGMENT_SLOTS {
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
}
