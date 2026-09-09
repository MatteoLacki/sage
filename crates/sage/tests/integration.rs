//! Ensure that we exhaustively visit all fragment ions matching tolerances

use quickcheck_macros::quickcheck;
use sage_core::database::{Builder, IndexedDatabase, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::mass::{Tolerance, PROTON};
use sage_core::scoring::{RankingScore, ScoreType, Scorer};
use sage_core::spectrum::{Peak, Precursor, ProcessedSpectrum};
use sage_core::spline::{Extrapolation, LinearSpline, ValueTolSpline};

const FASTA: &'static str = r#"
>sp|Q99536|VAT1_HUMAN Synaptic vesicle membrane protein VAT-1 homolog OS=Homo sapiens OX=9606 GN=VAT1 PE=1 SV=2
MSDEREVAEAATGEDASSPPPKTEAASDPQHPAASEGAAAAAASPPLLRCLVLTGFGGYD
KVKLQSRPAAPPAPGPGQLTLRLRACGLNFADLMARQGLYDRLPPLPVTPGMEGAGVVIA
VGEGVSDRKAGDRVMVLNRSGMWQEEVTVPSVQTFLIPEAMTFEEAAALLVNYITAYMVL
FDFGNLQPGHSVLVHMAAGGVGMAAVQLCRTVENVTVFGTASASKHEALKENGVTHPIDY
HTTDYVDEIKKISPKGVDIVMDPLGGSDTAKGYNLLKPMGKVVTYGMANLLTGPKRNLMA
LARTWWNQFSVTALQLLQANRAVCGFHLGYLDGEVELVSGVVARLLALYNQGHIKPHIDS
VWPFEKVADAMKQMQEKKNVGKVLLVPGPEKEN
"#;

fn mk_database(bucket_size: usize) -> IndexedDatabase {
    let builder = Builder {
        bucket_size: Some(bucket_size),
        fasta: Some("static".into()),
        ..Default::default()
    };
    let fasta = Fasta::parse(FASTA.into(), "rev_", false);

    builder.make_parameters().build(fasta)
}

#[quickcheck]
fn check_all_ions_visited(target_fragment_mz: f32, bucket_size: usize) {
    let database = mk_database(bucket_size.clamp(1, 8192));

    // Map PeptideIx -> number of fragments between 500 & 700 m/z
    // We want to make sure that IndexedDatabase::query hits *all* of them
    let mut expected = vec![0usize; database.peptides.len()];

    let fragment_tol = Tolerance::Da(-100.0, 100.0);
    let (frag_lo, frag_hi) = fragment_tol.bounds(target_fragment_mz);

    for (chunk_idx, chunk) in database.fragments.chunks(database.bucket_size).enumerate() {
        // Check for total ordering by PeptideIx within a chunk
        let mut last = PeptideIx(0);
        for frag in chunk {
            assert!(frag.peptide_index >= last);
            assert!(frag.fragment_mz >= database.buckets()[chunk_idx]);
            if chunk_idx + 1 < database.buckets().len() {
                assert!(frag.fragment_mz <= database.buckets()[chunk_idx + 1]);
            }

            if frag.fragment_mz >= frag_lo && frag.fragment_mz <= frag_hi {
                expected[frag.peptide_index.0 as usize] += 1;
            }
            last = frag.peptide_index;
        }
    }

    let mut visited = vec![0usize; database.peptides.len()];

    // Hit all peptides in database, track how many of the 500-700 fragment m/z's
    // are returned to us by searching the database.
    let query = database.query(1000.0, Tolerance::Da(-5000.0, 5000.0));

    for fragment in query.page_search(target_fragment_mz, fragment_tol) {
        visited[fragment.peptide_index.0 as usize] += 1;
    }

    // Make sure we visited every possible fragment
    assert_eq!(expected, visited);
}

/// Fasta for the per-precursor ppm tolerance tests: three clean tryptic
/// peptides of increasing mass (~600, ~945 Da), none is index 0 by mass so
/// filtering on `PeptideIx::default()` in `build_features` can't hide a
/// genuine hit. `PEPTIDEK` is deliberately placed first in the protein so
/// its cleavage site isn't blocked by the "no cleavage before P" rule.
const PPM_WINDOW_FASTA: &'static str = ">tol_test\nPEPTIDEKGGGGGGGGKMAAAAAAK\n";
const TARGET_SEQUENCE: &[u8] = b"PEPTIDEK";

fn mk_ppm_window_database() -> IndexedDatabase {
    let builder = Builder {
        bucket_size: Some(64),
        fasta: Some("static".into()),
        generate_decoys: Some(false),
        ..Default::default()
    };
    let fasta = Fasta::parse(PPM_WINDOW_FASTA.into(), "rev_", false);
    builder.make_parameters().build(fasta)
}

/// Real b/y fragment peaks for `TARGET_SEQUENCE`, taken directly from the
/// database's own theoretical fragments so they are guaranteed to match
/// under a tight fragment tolerance.
fn target_peaks(db: &IndexedDatabase) -> (PeptideIx, Vec<Peak>) {
    let target_idx = db
        .peptides
        .iter()
        .position(|p| &*p.sequence == TARGET_SEQUENCE)
        .expect("target peptide present in digest");
    let target_idx = PeptideIx(target_idx as u32);

    let peaks = db
        .fragments
        .iter()
        .filter(|f| f.peptide_index == target_idx)
        .map(|f| Peak {
            mass: f.fragment_mz,
            intensity: 100.0,
        })
        .collect();

    (target_idx, peaks)
}

fn mk_scorer(db: &IndexedDatabase, precursor_tol: Tolerance) -> Scorer<'_> {
    Scorer {
        db,
        precursor_tol,
        fragment_tol: Tolerance::Da(-0.01, 0.01),
        min_matched_peaks: 1,
        min_isotope_err: 0,
        max_isotope_err: 0,
        min_precursor_charge: 1,
        max_precursor_charge: 1,
        override_precursor_charge: false,
        max_fragment_charge: None,
        chimera: false,
        report_psms: 5,
        wide_window: false,
        predicted_rt: None,
        predicted_iim: None,
        rt_tol: None,
        mobility_tol: None,
        rt_sigma: None,
        iim_sigma: None,
        annotate_matches: false,
        score_type: ScoreType::SageHyperScore,
        ranking_score: RankingScore::CombinedScore,
        predicted_fragment_intensity_index: None,
        predicted_fragment_intensity_annotation_id: None,
        predicted_fragment_intensity: None,
    }
}

/// Build a `[M+H]+` m/z that is `offset_ppm` away from `monoisotopic`.
fn shifted_precursor_mz(monoisotopic: f32, offset_ppm: f32) -> f32 {
    let shifted_mass = monoisotopic * (1.0 + offset_ppm / 1_000_000.0);
    shifted_mass + PROTON
}

/// A candidate whose recorded precursor m/z is 30 ppm off from its
/// theoretical mass is unreachable under a narrow global `precursor_tol`
/// (±5 ppm) when the precursor carries no per-precursor override — this is
/// the pre-existing, unchanged behavior.
#[test]
fn candidate_unreachable_without_custom_window() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 30.0);

    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None, // no per-precursor override
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "no-window".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let scorer = mk_scorer(&db, Tolerance::Ppm(-5.0, 5.0));
    let features = scorer.score_standard(&query);
    assert!(
        features.is_empty(),
        "narrow global tol should not reach a candidate 30 ppm away"
    );
}

/// The same candidate as above IS reachable once the precursor carries its
/// own wide `isolation_window`, even though the run-global `precursor_tol`
/// is unchanged (and still too narrow on its own) — this exercises the
/// per-precursor ppm override added on top of `Precursor::isolation_window`.
#[test]
fn candidate_reachable_with_wide_custom_window() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 30.0);

    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: Some(Tolerance::Ppm(-50.0, 50.0)), // custom per-precursor override
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "wide-window".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    // Same narrow global tol as the "unreachable" test above.
    let scorer = mk_scorer(&db, Tolerance::Ppm(-5.0, 5.0));
    let features = scorer.score_standard(&query);
    assert_eq!(
        features.len(),
        1,
        "wide per-precursor window should reach the 30 ppm-off candidate"
    );
    assert_eq!(features[0].peptide_idx, target_idx);
}

/// A flat (value-independent) `ValueTolSpline`, `[lo, hi]` everywhere --
/// the "robust flat window" shape actually used for `rt_tol_sec`/
/// `mobility_tol` (a 2-node spline with identical values at both nodes,
/// not a separate flat-tolerance type).
fn flat_value_tol_spline(lo: f32, hi: f32) -> ValueTolSpline {
    let flat = |value: f32| LinearSpline {
        grid_start: 0.0,
        grid_step: 1000.0,
        values: vec![value, value],
        extrapolation: Extrapolation::Flat,
    };
    ValueTolSpline {
        lo: flat(lo),
        hi: flat(hi),
    }
}

/// A flat (value-independent) `LinearSpline` -- same role for `rt_sigma`/
/// `iim_sigma`-style scale fields as `flat_value_tol_spline` has for
/// `rt_tol`/`mobility_tol`, now that `rt_sigma` is spline-shaped too
/// (`plans/rt_heteroscedastic_tolerance_spline.md`).
fn flat_linear_spline(value: f32) -> LinearSpline {
    LinearSpline {
        grid_start: 0.0,
        grid_step: 1000.0,
        values: vec![value, value],
        extrapolation: Extrapolation::Flat,
    }
}

/// A candidate whose entry in `predicted_rt` gives it a predicted RT far
/// outside `rt_tol` of the observed spectrum's `scan_start_time` is
/// unreachable, even though precursor/fragment tolerances alone would have
/// matched it fine (same query as the RT-mismatch test below, minus the
/// `predicted_rt` map, does find it — see that test). `predicted_iim` is
/// left unconfigured entirely, confirming RT filtering works independently
/// of IIM — see plans/rt_iim_independent_dimensions.md.
#[test]
fn candidate_unreachable_outside_rt_tol() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "rt-mismatch".into(),
        scan_start_time: 10.0,
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_rt = vec![None; db.peptides.len()];
    predicted_rt[target_idx.0 as usize] = Some(15.0f32); // far outside [9.8, 10.2]

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_rt = Some(&predicted_rt);
    scorer.rt_tol = Some(flat_value_tol_spline(-0.2, 0.2));
    scorer.predicted_iim = None;
    scorer.mobility_tol = None;

    let features = scorer.score_standard(&query);
    assert!(
        features.is_empty(),
        "predicted RT 5 minutes outside rt_tol should reject the candidate"
    );
}

/// Same query, target, and `predicted_rt` map shape as above, but with a
/// predicted RT inside `rt_tol` — confirms the candidate is reachable (i.e.
/// the rejection above is really coming from the RT mismatch, not some
/// other difference), and that no `predicted_iim`/`mobility_tol` at all
/// (not even a `None` observed `inverse_ion_mobility`, the dimension is
/// simply not configured) never touches eviction.
#[test]
fn candidate_reachable_within_rt_tol() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "rt-match".into(),
        scan_start_time: 10.0,
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_rt = vec![None; db.peptides.len()];
    predicted_rt[target_idx.0 as usize] = Some(10.05f32); // inside [9.8, 10.2]

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_rt = Some(&predicted_rt);
    scorer.rt_tol = Some(flat_value_tol_spline(-0.2, 0.2));
    scorer.predicted_iim = None;
    scorer.mobility_tol = None;

    let features = scorer.score_standard(&query);
    assert_eq!(
        features.len(),
        1,
        "predicted RT inside rt_tol should reach the candidate"
    );
    assert_eq!(features[0].peptide_idx, target_idx);
}

/// External-prediction z² feature (`plans/lda_external_rt_iim_features.md`):
/// same reachable-candidate setup as `candidate_reachable_within_rt_tol`,
/// plus `rt_sigma` — confirms `Feature::predicted_rt_external`/
/// `delta_rt_z2_external` are populated from the same dense array
/// `evict_rt_iim_mismatches` reads, independent of SAGE's own internal
/// `predicted_rt`/`delta_rt_model` (which `score_standard` alone, without
/// the `retention_model`/`mobility_model` post-processing pass runner.rs
/// applies, never touches — left at their `Default` `0.0`/stub here).
#[test]
fn delta_rt_z2_external_computed_from_sigma() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "rt-z2".into(),
        scan_start_time: 10.0,
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_rt = vec![None; db.peptides.len()];
    predicted_rt[target_idx.0 as usize] = Some(10.05f32);

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_rt = Some(&predicted_rt);
    scorer.rt_tol = Some(flat_value_tol_spline(-0.2, 0.2));
    scorer.rt_sigma = Some(flat_linear_spline(0.1)); // minutes, same unit as rt_tol's converted values

    let features = scorer.score_standard(&query);
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].predicted_rt_external, 10.05);
    // z = (10.0 - 10.05) / 0.1 = -0.5, z^2 = 0.25
    assert!(
        (features[0].delta_rt_z2_external - 0.25).abs() < 1e-5,
        "expected z^2 ~= 0.25, got {}",
        features[0].delta_rt_z2_external
    );
}

/// `rt_sigma` left unset even though `predicted_rt`/`rt_tol` are configured
/// (shouldn't happen via `Input::build`'s validation, but `Scorer` itself
/// doesn't enforce the pairing) — the external z² feature stays at its
/// `Default` `0.0` rather than dividing by a missing/zero sigma.
#[test]
fn delta_rt_z2_external_zero_without_sigma() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "rt-z2-no-sigma".into(),
        scan_start_time: 10.0,
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_rt = vec![None; db.peptides.len()];
    predicted_rt[target_idx.0 as usize] = Some(10.05f32);

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_rt = Some(&predicted_rt);
    scorer.rt_tol = Some(flat_value_tol_spline(-0.2, 0.2));
    // scorer.rt_sigma left None

    let features = scorer.score_standard(&query);
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].predicted_rt_external, 0.0);
    assert_eq!(features[0].delta_rt_z2_external, 0.0);
}

/// IIM-only filtering, the symmetric counterpart of the RT tests above —
/// `predicted_rt`/`rt_tol` are left entirely unconfigured, confirming IIM
/// filtering also works independently.
#[test]
fn candidate_unreachable_outside_iim_tol() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: Some(1.0),
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "iim-mismatch".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_iim = vec![None; db.peptides.len()];
    let slot = Scorer::iim_dense_slot(target_idx.0 as usize, 1, 1, 1).unwrap();
    predicted_iim[slot] = Some(1.5f32); // far outside [0.9, 1.1]

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_iim = Some(&predicted_iim);
    scorer.mobility_tol = Some(flat_value_tol_spline(-0.1, 0.1));
    scorer.predicted_rt = None;
    scorer.rt_tol = None;

    let features = scorer.score_standard(&query);
    assert!(
        features.is_empty(),
        "predicted IIM far outside mobility_tol should reject the candidate"
    );
}

/// Same as above but with a predicted IIM inside `mobility_tol`.
#[test]
fn candidate_reachable_within_iim_tol() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: Some(1.0),
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "iim-match".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_iim = vec![None; db.peptides.len()];
    let slot = Scorer::iim_dense_slot(target_idx.0 as usize, 1, 1, 1).unwrap();
    predicted_iim[slot] = Some(1.05f32); // inside [0.9, 1.1]

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_iim = Some(&predicted_iim);
    scorer.mobility_tol = Some(flat_value_tol_spline(-0.1, 0.1));
    scorer.predicted_rt = None;
    scorer.rt_tol = None;

    let features = scorer.score_standard(&query);
    assert_eq!(
        features.len(),
        1,
        "predicted IIM inside mobility_tol should reach the candidate"
    );
    assert_eq!(features[0].peptide_idx, target_idx);
}

/// IIM counterpart of `delta_rt_z2_external_computed_from_sigma` — same
/// shape, `Scorer::iim_dense_slot` indexing instead of a flat per-peptide
/// array.
#[test]
fn delta_ims_z2_external_computed_from_sigma() {
    let db = mk_ppm_window_database();
    let (target_idx, peaks) = target_peaks(&db);
    assert!(!peaks.is_empty(), "expected real fragment peaks for target");

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: Some(1.0),
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "iim-z2".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_iim = vec![None; db.peptides.len()];
    let slot = Scorer::iim_dense_slot(target_idx.0 as usize, 1, 1, 1).unwrap();
    predicted_iim[slot] = Some(1.05f32);

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.predicted_iim = Some(&predicted_iim);
    scorer.mobility_tol = Some(flat_value_tol_spline(-0.1, 0.1));
    scorer.iim_sigma = Some(0.02); // unitless (1/K0), no conversion

    let features = scorer.score_standard(&query);
    assert_eq!(features.len(), 1);
    assert_eq!(features[0].predicted_ims_external, 1.05);
    // z = (1.0 - 1.05) / 0.02 = -2.5, z^2 = 6.25
    assert!(
        (features[0].delta_ims_z2_external - 6.25).abs() < 1e-4,
        "expected z^2 ~= 6.25, got {}",
        features[0].delta_ims_z2_external
    );
}

/// Two proteins whose sole tryptic peptides differ only by an I<->L swap —
/// isomeric (identical monoisotopic residue mass), so both have byte-identical
/// theoretical fragment masses and therefore tie exactly on hyperscore for any
/// shared peak list. This isolates the new `combined_score` ranking (`build_features`
/// in `scoring.rs`) from hyperscore itself: with a real hyperscore tie, RT
/// closeness is the only thing that can break it.
const ISOBARIC_FASTA: &'static str = ">iso1\nMPEPTIDEK\n>iso2\nMPEPTLDEK\n";

fn mk_isobaric_database() -> IndexedDatabase {
    let builder = Builder {
        bucket_size: Some(64),
        fasta: Some("static".into()),
        generate_decoys: Some(false),
        ..Default::default()
    };
    let fasta = Fasta::parse(ISOBARIC_FASTA.into(), "rev_", false);
    builder.make_parameters().build(fasta)
}

/// A candidate closer to its predicted RT outranks a hyperscore-tied
/// candidate that's farther from its predicted RT — the soft `combined_score`
/// penalty (`hyperscore - 0.5*(z_rt_external^2 + z_iim_external^2)`) reorders
/// otherwise-equal candidates instead of hard-evicting either one.
#[test]
fn combined_score_ranks_rt_tied_hyperscore_candidates() {
    let db = mk_isobaric_database();
    let iso1 = PeptideIx(
        db.peptides
            .iter()
            .position(|p| &*p.sequence == b"MPEPTIDEK")
            .expect("iso1 present in digest") as u32,
    );
    let iso2 = PeptideIx(
        db.peptides
            .iter()
            .position(|p| &*p.sequence == b"MPEPTLDEK")
            .expect("iso2 present in digest") as u32,
    );
    assert_eq!(
        db[iso1].monoisotopic, db[iso2].monoisotopic,
        "I<->L substitution must be exactly isomeric for this test to isolate RT"
    );

    let peaks: Vec<Peak> = db
        .fragments
        .iter()
        .filter(|f| f.peptide_index == iso1)
        .map(|f| Peak {
            mass: f.fragment_mz,
            intensity: 100.0,
        })
        .collect();
    assert!(!peaks.is_empty(), "expected real fragment peaks for iso1");

    let mz = shifted_precursor_mz(db[iso1].monoisotopic, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "rt-tiebreak".into(),
        scan_start_time: 10.0,
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    // iso1 predicted far from observed RT (z^2 large), iso2 predicted right
    // on top of it (z^2 ~= 0) -- both well inside the wide rt_tol, so neither
    // is evicted; only the ranking should differ.
    let mut predicted_rt = vec![None; db.peptides.len()];
    predicted_rt[iso1.0 as usize] = Some(11.0f32);
    predicted_rt[iso2.0 as usize] = Some(10.0f32);

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.report_psms = 2;
    scorer.predicted_rt = Some(&predicted_rt);
    scorer.rt_tol = Some(flat_value_tol_spline(-5.0, 5.0));
    scorer.rt_sigma = Some(flat_linear_spline(0.1));

    let features = scorer.score_standard(&query);
    assert_eq!(features.len(), 2, "both isomeric candidates should survive");
    assert_eq!(
        features[0].hyperscore, features[1].hyperscore,
        "isomeric candidates must tie exactly on raw hyperscore"
    );
    assert_eq!(
        features[0].peptide_idx, iso2,
        "candidate closer to predicted RT should rank first despite the hyperscore tie"
    );
    // z_iso1 = (10.0 - 11.0) / 0.1 = -10, term = -0.5*100 = -50; z_iso2 = 0.
    // delta_next = combined_best - combined_next = H - (H - 50) = 50, not the
    // 0.0 a straight hyperscore-tie margin would give.
    assert!(
        (features[0].delta_next - 50.0).abs() < 1e-3,
        "delta_next should reflect the combined_score gap, not a 0.0 hyperscore tie: got {}",
        features[0].delta_next
    );
}

/// Same isomeric-tie setup as `combined_score_ranks_rt_tied_hyperscore_candidates`,
/// but with `ranking_score: RankingScore::Hyperscore` -- the runtime-selectable
/// escape hatch back to pre-`combined_score` behavior (config `ranking_score`,
/// same `Option<T>` + `Input::build()`-default shape as `score_type`). With a
/// real hyperscore tie and no RT influence on ranking, `build_features` can't
/// break the tie: `delta_next` stays exactly `0.0` regardless of how far apart
/// the two candidates' predicted RTs are.
#[test]
fn ranking_score_hyperscore_ignores_rt_penalty() {
    let db = mk_isobaric_database();
    let iso1 = PeptideIx(
        db.peptides
            .iter()
            .position(|p| &*p.sequence == b"MPEPTIDEK")
            .expect("iso1 present in digest") as u32,
    );
    let iso2 = PeptideIx(
        db.peptides
            .iter()
            .position(|p| &*p.sequence == b"MPEPTLDEK")
            .expect("iso2 present in digest") as u32,
    );

    let peaks: Vec<Peak> = db
        .fragments
        .iter()
        .filter(|f| f.peptide_index == iso1)
        .map(|f| Peak {
            mass: f.fragment_mz,
            intensity: 100.0,
        })
        .collect();
    assert!(!peaks.is_empty(), "expected real fragment peaks for iso1");

    let mz = shifted_precursor_mz(db[iso1].monoisotopic, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        inverse_ion_mobility: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "rt-tiebreak-hyperscore-mode".into(),
        scan_start_time: 10.0,
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let mut predicted_rt = vec![None; db.peptides.len()];
    predicted_rt[iso1.0 as usize] = Some(11.0f32);
    predicted_rt[iso2.0 as usize] = Some(10.0f32);

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0));
    scorer.report_psms = 2;
    scorer.predicted_rt = Some(&predicted_rt);
    scorer.rt_tol = Some(flat_value_tol_spline(-5.0, 5.0));
    scorer.rt_sigma = Some(flat_linear_spline(0.1));
    scorer.ranking_score = RankingScore::Hyperscore;

    let features = scorer.score_standard(&query);
    assert_eq!(features.len(), 2, "both isomeric candidates should survive");
    assert_eq!(
        features[0].hyperscore, features[1].hyperscore,
        "isomeric candidates must tie exactly on raw hyperscore"
    );
    assert_eq!(
        features[0].delta_next, 0.0,
        "Hyperscore mode must ignore the RT penalty entirely -- a real tie stays a tie"
    );
}
