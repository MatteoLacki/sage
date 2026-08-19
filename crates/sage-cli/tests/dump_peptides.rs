use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use sage_core::database::Builder;
use std::collections::HashSet;
use std::fs::File;
use std::process::Command;

/// Cross-checks the standalone `dump_peptides` binary against calling the
/// digestion library directly with the same `Builder::default()` settings —
/// confirms the binary is a faithful, fasta-only proxy for
/// `Parameters::digest` (see `plans/better_sage_filtering.md`), not just
/// that it runs.
#[test]
fn dump_peptides_matches_direct_digestion() -> anyhow::Result<()> {
    let fasta_path = "../../tests/trivial_3proteins.fasta";

    // Expected: call the digestion library directly, same defaults the
    // binary uses when no `--config` is given.
    let mut builder = Builder::default();
    builder.update_fasta(fasta_path.into());
    let fasta = sage_cloudpath::util::read_fasta(fasta_path, "rev_", true)?;
    let expected_peptides = builder.make_parameters().digest(&fasta);
    assert!(
        !expected_peptides.is_empty(),
        "fixture should digest to at least one peptide"
    );
    let expected: HashSet<(String, bool)> = expected_peptides
        .iter()
        .map(|p| (p.to_string(), p.decoy))
        .collect();
    assert!(
        expected.iter().any(|(_, decoy)| !decoy),
        "expected at least one target peptide"
    );
    assert!(
        expected.iter().any(|(_, decoy)| *decoy),
        "expected at least one decoy peptide"
    );

    // Actual: invoke the compiled binary.
    let output_path = std::env::temp_dir().join(format!(
        "dump_peptides_test_{}_peptides.parquet",
        std::process::id()
    ));

    let status = Command::new(env!("CARGO_BIN_EXE_dump_peptides"))
        .arg("-f")
        .arg(fasta_path)
        .arg("-o")
        .arg(&output_path)
        .status()?;
    assert!(status.success(), "dump_peptides exited non-zero");

    let file = File::open(&output_path)?;
    let reader = SerializedFileReader::new(file)?;
    let schema = reader.metadata().file_metadata().schema_descr();
    let col_idx = |name: &str| schema.columns().iter().position(|c| c.name() == name).unwrap();
    let (i_peptide, i_decoy) = (col_idx("peptide"), col_idx("decoy"));

    let actual: HashSet<(String, bool)> = reader
        .get_row_iter(None)?
        .map(|row| {
            let row = row.unwrap();
            (row.get_string(i_peptide).unwrap().clone(), row.get_bool(i_decoy).unwrap())
        })
        .collect();

    let _ = std::fs::remove_file(&output_path);

    assert_eq!(
        actual, expected,
        "dump_peptides output diverges from direct Parameters::digest"
    );

    Ok(())
}
