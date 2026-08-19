//! Reader for externally-predicted peptide RT/IIM properties (see
//! `necromerge2`'s `git/featureprediction` and `plans/better_sage_filtering.md`).
//!
//! One parquet file, columns `sequence`(utf8, `Peptide::Display` format),
//! `charge`(int32), `rt`(double, minutes), `iim`(double, 1/K0) — one row per
//! `(sequence, charge)` pair. Loaded once at startup into a plain
//! `HashMap<(String, u8), (f32, f32)>`.

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

pub fn read_predicted_properties(
    path: &Path,
) -> Result<HashMap<(String, u8), (f32, f32)>, PredictedPropertiesError> {
    let file = File::open(path)?;
    let reader = SerializedFileReader::new(file)?;

    let schema = reader.metadata().file_metadata().schema_descr();
    let col_idx = |name: &'static str| -> Result<usize, PredictedPropertiesError> {
        schema
            .columns()
            .iter()
            .position(|c| c.name() == name)
            .ok_or(PredictedPropertiesError::MissingColumn(name))
    };

    let i_sequence = col_idx("sequence")?;
    let i_charge = col_idx("charge")?;
    let i_rt = col_idx("rt")?;
    let i_iim = col_idx("iim")?;

    let mut map = HashMap::new();
    for row in reader.get_row_iter(None)? {
        let row = row?;
        let sequence = row.get_string(i_sequence)?.clone();
        let charge = row.get_int(i_charge)? as u8;
        let rt = row.get_double(i_rt)? as f32;
        let iim = row.get_double(i_iim)? as f32;
        map.insert((sequence, charge), (rt, iim));
    }

    Ok(map)
}
