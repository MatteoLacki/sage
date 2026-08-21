//! Readers for externally-predicted peptide RT/IIM properties (see
//! `necromerge2`'s `git/featureprediction` and
//! `plans/rt_iim_independent_dimensions.md`). RT and IIM are two
//! independent, separately-loadable files -- RT has no charge dimension
//! (Chronologer_RT doesn't take one), IIM does.

#![cfg(feature = "parquet")]

use std::collections::HashMap;
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

/// One parquet file, columns `sequence`(utf8, `Peptide::Display` format),
/// `rt`(double, minutes) -- one row per sequence.
pub fn read_predicted_rt(path: &Path) -> Result<HashMap<String, f32>, PredictedPropertiesError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    let schema = reader.metadata().file_metadata().schema_descr();

    let i_sequence = col_idx(schema, "sequence")?;
    let i_rt = col_idx(schema, "rt")?;

    let mut map = HashMap::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        let sequence = row.get_string(i_sequence)?.clone();
        let rt = row.get_double(i_rt)? as f32;
        map.insert(sequence, rt);
    }

    Ok(map)
}

/// One parquet file, columns `sequence`(utf8, `Peptide::Display` format),
/// `charge`(int32), `iim`(double, 1/K0) -- one row per `(sequence, charge)`
/// pair.
pub fn read_predicted_iim(
    path: &Path,
) -> Result<HashMap<(String, u8), f32>, PredictedPropertiesError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;
    let schema = reader.metadata().file_metadata().schema_descr();

    let i_sequence = col_idx(schema, "sequence")?;
    let i_charge = col_idx(schema, "charge")?;
    let i_iim = col_idx(schema, "iim")?;

    let mut map = HashMap::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        let sequence = row.get_string(i_sequence)?.clone();
        let charge = row.get_int(i_charge)? as u8;
        let iim = row.get_double(i_iim)? as f32;
        map.insert((sequence, charge), iim);
    }

    Ok(map)
}
