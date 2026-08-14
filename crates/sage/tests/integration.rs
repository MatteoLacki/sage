//! Ensure that we exhaustively visit all fragment ions matching tolerances

use quickcheck_macros::quickcheck;
use sage_core::database::{Builder, IndexedDatabase, PeptideIx};
use sage_core::fasta::Fasta;
use sage_core::mass::{Tolerance, PROTON};
use sage_core::scoring::{ScoreType, Scorer};
use sage_core::spectrum::{Peak, Precursor, ProcessedSpectrum};

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
    let query = database.query(1000.0, Tolerance::Da(-5000.0, 5000.0), fragment_tol);

    for fragment in query.page_search(target_fragment_mz) {
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
