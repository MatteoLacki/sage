//! Readers for externally-predicted peptide RT/IIM properties (see
//! `necromerge2`'s `git/featureprediction` and
//! `plans/rt_iim_independent_dimensions.md`). RT and IIM are two
//! independent, separately-loadable files -- RT has no charge dimension
//! (Chronologer_RT doesn't take one), IIM does.
//!
//! Both files are positional, not `sequence`-keyed (2026-09-08 -- see
//! `docs/ai/dumped_peptides_positional_predictions.md`): row `i` of
//! `predicted_rt.parquet` is peptide row `i` of the `dumped_peptides`
//! parquet both files were derived from; `predicted_iim.parquet` is dense
//! over `(peptide_row, charge)` in the same `peptide_idx * charge_span +
//! (charge - min_charge)` layout `Scorer::iim_dense_slot` already uses.
//! Each file carries a `dumped_peptides_sha256` file-level key-value
//! metadata entry (see `dumped_peptides_fingerprint`) -- the caller
//! (`sage-cli`'s `resolve_predicted_rt`/`resolve_predicted_iim`) must check
//! it against the current run's own digest before trusting these
//! positions.

use std::fs::File;
use std::path::Path;

use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;

#[derive(Debug, thiserror::Error)]
pub enum PredictedPropertiesError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parquet error: {0}")]
    Parquet(#[from] parquet::errors::ParquetError),
    #[error("missing column: {0}")]
    MissingColumn(&'static str),
    #[error("{0} is missing required `dumped_peptides_sha256` file metadata -- regenerate it with a current `feature-prediction` build")]
    MissingFingerprint(String),
}

fn col_idx(
    schema: &parquet::schema::types::SchemaDescriptor,
    name: &'static str,
) -> Result<usize, PredictedPropertiesError> {
    schema
        .columns()
        .iter()
        .position(|c| c.name() == name)
        .ok_or(PredictedPropertiesError::MissingColumn(name))
}

fn read_fingerprint(
    reader: &SerializedFileReader<File>,
    path: &Path,
) -> Result<String, PredictedPropertiesError> {
    reader
        .metadata()
        .file_metadata()
        .key_value_metadata()
        .into_iter()
        .flatten()
        .find(|kv| kv.key == "dumped_peptides_sha256")
        .and_then(|kv| kv.value.clone())
        .ok_or_else(|| PredictedPropertiesError::MissingFingerprint(path.display().to_string()))
}

/// One parquet file, single column `rt`(double, minutes), one row per
/// `dumped_peptides` row, in that same order. Returns the embedded
/// `dumped_peptides` fingerprint alongside the dense values.
pub fn read_predicted_rt(
    path: &Path,
) -> Result<(String, Vec<f32>), PredictedPropertiesError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    let fingerprint = read_fingerprint(&reader, path)?;
    let schema = reader.metadata().file_metadata().schema_descr();
    let i_rt = col_idx(schema, "rt")?;

    let mut values = Vec::with_capacity(reader.metadata().file_metadata().num_rows() as usize);
    for row in reader.get_row_iter(None)? {
        let row = row?;
        values.push(row.get_double(i_rt)? as f32);
    }

    Ok((fingerprint, values))
}

/// One parquet file, single column `iim`(double, 1/K0), dense over
/// `(peptide_row, charge)` -- see module docs for the slot layout. Returns
/// the embedded `dumped_peptides` fingerprint alongside the dense values.
pub fn read_predicted_iim(
    path: &Path,
) -> Result<(String, Vec<f32>), PredictedPropertiesError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    let fingerprint = read_fingerprint(&reader, path)?;
    let schema = reader.metadata().file_metadata().schema_descr();
    let i_iim = col_idx(schema, "iim")?;

    let mut values = Vec::with_capacity(reader.metadata().file_metadata().num_rows() as usize);
    for row in reader.get_row_iter(None)? {
        let row = row?;
        values.push(row.get_double(i_iim)? as f32);
    }

    Ok((fingerprint, values))
}
