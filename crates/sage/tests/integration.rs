//! Ensure that we exhaustively visit all fragment ions matching tolerances

use quickcheck_macros::quickcheck;
use sage_core::database::{Builder, IndexedDatabase, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::mass::{Tolerance, PROTON};
use sage_core::scoring::{ScoreType, Scorer};
use sage_core::spectrum::{Peak, Precursor, ProcessedSpectrum};
use sage_core::spline::{FragmentTolSpline, LinearSpline};

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

fn mk_scorer(
    db: &IndexedDatabase,
    precursor_tol: Tolerance,
    fragment_tol_spline: Option<FragmentTolSpline>,
) -> Scorer<'_> {
    mk_scorer_with_fragment_tol(db, precursor_tol, Tolerance::Da(-0.01, 0.01), fragment_tol_spline)
}

fn mk_scorer_with_fragment_tol(
    db: &IndexedDatabase,
    precursor_tol: Tolerance,
    fragment_tol: Tolerance,
    fragment_tol_spline: Option<FragmentTolSpline>,
) -> Scorer<'_> {
    Scorer {
        db,
        precursor_tol,
        fragment_tol,
        fragment_tol_spline,
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
        predicted_properties: None,
        rt_tol: None,
        mobility_tol: None,
        annotate_matches: false,
        score_type: ScoreType::SageHyperScore,
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

    let scorer = mk_scorer(&db, Tolerance::Ppm(-5.0, 5.0), None);
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
    let scorer = mk_scorer(&db, Tolerance::Ppm(-5.0, 5.0), None);
    let features = scorer.score_standard(&query);
    assert_eq!(
        features.len(),
        1,
        "wide per-precursor window should reach the 30 ppm-off candidate"
    );
    assert_eq!(features[0].peptide_idx, target_idx);
}

/// Shift every fragment peak's mass by `offset_ppm`, simulating a spectrum
/// whose fragments are all systematically off by a constant relative error
/// (what a mass-dependent calibration spline is meant to correct for).
fn shift_peaks(peaks: &[Peak], offset_ppm: f32) -> Vec<Peak> {
    peaks
        .iter()
        .map(|p| Peak {
            mass: p.mass * (1.0 + offset_ppm / 1_000_000.0),
            intensity: p.intensity,
        })
        .collect()
}

/// A flat (constant-valued) spline covering any plausible fragment mass in
/// these tests, ±`ppm` on both edges.
fn flat_fragment_tol_spline(ppm: f32) -> FragmentTolSpline {
    let flat = |sign: f32| LinearSpline {
        grid_start: 0.0,
        grid_step: 2000.0,
        values: vec![sign * ppm, sign * ppm],
    };
    FragmentTolSpline {
        ppm_lo: flat(-1.0),
        ppm_hi: flat(1.0),
    }
}

/// All fragment peaks shifted 100 ppm off their theoretical masses are
/// unreachable under a narrow flat `fragment_tol` (±5 ppm) with no spline —
/// this is the pre-existing, unchanged behavior.
#[test]
fn candidate_unreachable_without_fragment_tol_spline() {
    let db = mk_ppm_window_database();
    let (target_idx, real_peaks) = target_peaks(&db);
    assert!(
        !real_peaks.is_empty(),
        "expected real fragment peaks for target"
    );
    let peaks = shift_peaks(&real_peaks, 100.0);

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "no-spline".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    let scorer = mk_scorer_with_fragment_tol(
        &db,
        Tolerance::Ppm(-50.0, 50.0),
        Tolerance::Ppm(-5.0, 5.0),
        None,
    );
    let features = scorer.score_standard(&query);
    assert!(
        features.is_empty(),
        "narrow flat fragment_tol should not reach fragments 100 ppm away"
    );
}

/// The same 100 ppm-shifted fragments ARE reachable once a
/// `fragment_tol_spline` wide enough to cover the shift is configured, even
/// though the flat `fragment_tol` is unchanged (and still too narrow on its
/// own) — this exercises `Scorer::fragment_tol_spline` overriding the flat
/// tolerance per observed peak.
#[test]
fn candidate_reachable_with_fragment_tol_spline() {
    let db = mk_ppm_window_database();
    let (target_idx, real_peaks) = target_peaks(&db);
    assert!(
        !real_peaks.is_empty(),
        "expected real fragment peaks for target"
    );
    let peaks = shift_peaks(&real_peaks, 100.0);

    let target_mass = db[target_idx].monoisotopic;
    let mz = shifted_precursor_mz(target_mass, 0.0);
    let precursor = Precursor {
        mz,
        charge: Some(1),
        isolation_window: None,
        ..Default::default()
    };
    let query = ProcessedSpectrum {
        level: 2,
        id: "with-spline".into(),
        precursors: vec![precursor],
        peaks,
        ..Default::default()
    };

    // Same narrow flat fragment_tol as the "unreachable" test above.
    let scorer = mk_scorer_with_fragment_tol(
        &db,
        Tolerance::Ppm(-50.0, 50.0),
        Tolerance::Ppm(-5.0, 5.0),
        Some(flat_fragment_tol_spline(150.0)),
    );
    let features = scorer.score_standard(&query);
    assert_eq!(
        features.len(),
        1,
        "fragment_tol_spline covering the 100 ppm shift should reach the candidate"
    );
    assert_eq!(features[0].peptide_idx, target_idx);
}

/// A candidate whose entry in `predicted_properties` gives it a predicted RT
/// far outside `rt_tol` of the observed spectrum's `scan_start_time` is
/// unreachable, even though precursor/fragment tolerances alone would have
/// matched it fine (same query as the RT-mismatch test below, minus the
/// `predicted_properties` map, does find it — see that test).
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

    let key = (db[target_idx].to_string(), 1u8);
    let mut predicted = std::collections::HashMap::new();
    predicted.insert(key, (15.0f32, 0.0f32)); // rt=15.0, far outside [9.8, 10.2]

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0), None);
    scorer.predicted_properties = Some(&predicted);
    scorer.rt_tol = Some(Tolerance::Da(-0.2, 0.2));
    scorer.mobility_tol = None;

    let features = scorer.score_standard(&query);
    assert!(
        features.is_empty(),
        "predicted RT 5 minutes outside rt_tol should reject the candidate"
    );
}

/// Same query, target, and `predicted_properties` map shape as above, but
/// with a predicted RT inside `rt_tol` — confirms the candidate is reachable
/// (i.e. the rejection above is really coming from the RT mismatch, not
/// some other difference), and that a `None` observed `inverse_ion_mobility`
/// skips the IIM check entirely (no `mobility_tol` needed).
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

    let key = (db[target_idx].to_string(), 1u8);
    let mut predicted = std::collections::HashMap::new();
    predicted.insert(key, (10.05f32, 0.0f32)); // rt=10.05, inside [9.8, 10.2]

    let mut scorer = mk_scorer(&db, Tolerance::Ppm(-50.0, 50.0), None);
    scorer.predicted_properties = Some(&predicted);
    scorer.rt_tol = Some(Tolerance::Da(-0.2, 0.2));
    scorer.mobility_tol = None;

    let features = scorer.score_standard(&query);
    assert_eq!(
        features.len(),
        1,
        "predicted RT inside rt_tol should reach the candidate"
    );
    assert_eq!(features[0].peptide_idx, target_idx);
}
