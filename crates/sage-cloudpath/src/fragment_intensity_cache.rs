//! Reader for `git/featureprediction`'s fragment-intensity `PredictionCache`
//! (see `necromerge2`'s `docs/ai/predicted_fragment_intensity.md`).
//!
//! Two separate inputs, given as explicit paths:
//!   - a small, job-scoped pointer index, positional not `sequence`-keyed
//!     since 2026-09-08 (see `docs/ai/dumped_peptides_positional_predictions.md`):
//!     an mmappet dataset with `start`/`end` (int64, sentinel -1 for "no cache
//!     entry") -- half-open row range into the shared arrays below, dense over
//!     `(peptide_row, charge)` in the same layout `Scorer::iim_dense_slot`
//!     uses -- mapped in place. Produced by
//!     `feature-prediction-export-fragments-for-sage`. Its `meta.json` carries a
//!     `dumped_peptides_sha256` the caller must check (see
//!     `dumped_peptides_fingerprint`).
//!   - the shared, ever-growing `arrays.mmappet/` directory (schema.txt +
//!     `0.bin`/`1.bin`, same on-disk format `pmsms.rs` reads for
//!     `pmsms.mmappet`/precursors `.mmappet` -- schema/mmap-slice helpers
//!     are duplicated here, not shared with `pmsms.rs`, since that module's
//!     helpers are private and tied to its own `PmsmsError`; both are small,
//!     easily-tested, and deliberately kept apart per this project's
//!     documented-literal-duplication convention rather than adding a
//!     shared abstraction for two call sites).

use std::fs::File;
use std::path::Path;

use half::f16;
use memmap2::Mmap;

use crate::predicted_properties::{read_fingerprint, MappedColumn, PredictedPropertiesError};

#[derive(Debug, thiserror::Error)]
pub enum FragmentIntensityCacheError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error(transparent)]
    Predicted(#[from] PredictedPropertiesError),
    #[error("missing column: {0}")]
    MissingColumn(&'static str),
    #[error("mmappet column {0} has dtype `{1}`, expected `{2}`")]
    WrongDtype(&'static str, String, &'static str),
    #[error("{0}: `start` and `end` differ in length")]
    UnequalColumns(String),
}

/// The job-scoped, dense `(peptide_row, charge)` pointer index, mapped -- see
/// module docs. Values are permanent once written (the shared cache is
/// append-only).
#[derive(Clone)]
pub struct FragmentIntensityIndex {
    pub start: MappedColumn<i64>,
    pub end: MappedColumn<i64>,
}

/// Map the pointer index in `dir`. Returns the `dumped_peptides` fingerprint
/// alongside it -- the caller must check it against the current run's own digest
/// before trusting these positions.
pub fn read_fragment_intensity_index(
    dir: &Path,
) -> Result<(String, FragmentIntensityIndex), FragmentIntensityCacheError> {
    let fingerprint = read_fingerprint(dir)?;
    let start = MappedColumn::open(dir, "start")?;
    let end = MappedColumn::open(dir, "end")?;
    if start.len() != end.len() {
        return Err(FragmentIntensityCacheError::UnequalColumns(
            dir.display().to_string(),
        ));
    }
    Ok((fingerprint, FragmentIntensityIndex { start, end }))
}

/// Parse an mmappet `schema.txt` into an ordered list of `(dtype, column
/// name)`, matching the numpy dtype strings written by `git/mmappet`'s
/// Python writer.
fn read_mmappet_schema(dir: &Path) -> Result<Vec<(String, String)>, FragmentIntensityCacheError> {
    let text = std::fs::read_to_string(dir.join("schema.txt"))?;
    Ok(text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut parts = line.split_whitespace();
            let dtype = parts.next().unwrap_or_default().to_string();
            let name = parts.next().unwrap_or_default().to_string();
            (dtype, name)
        })
        .collect())
}

fn mmappet_column_index(
    schema: &[(String, String)],
    name: &'static str,
    expected_dtype: &'static str,
) -> Result<usize, FragmentIntensityCacheError> {
    let (idx, (dtype, _)) = schema
        .iter()
        .enumerate()
        .find(|(_, (_, col_name))| col_name == name)
        .ok_or(FragmentIntensityCacheError::MissingColumn(name))?;
    if dtype != expected_dtype {
        return Err(FragmentIntensityCacheError::WrongDtype(
            name,
            dtype.clone(),
            expected_dtype,
        ));
    }
    Ok(idx)
}

/// SAFETY: caller must have already checked `T`'s size/alignment matches the
/// on-disk dtype (verified via `mmappet_column_index`'s dtype check above);
/// the returned slice's lifetime is tied to the returned `Mmap`, which the
/// caller must keep alive for as long as the slice is used.
unsafe fn mmap_as_slice<T: Copy>(
    path: &Path,
) -> Result<(Mmap, &'static [T]), FragmentIntensityCacheError> {
    let file = File::open(path)?;
    let mmap = Mmap::map(&file)?;
    let n = mmap.len() / std::mem::size_of::<T>();
    let ptr = mmap.as_ptr() as *const T;
    let slice = std::slice::from_raw_parts(ptr, n);
    Ok((mmap, slice))
}

/// The shared, ever-growing `arrays.mmappet` dataset -- two flat mmapped
/// columns (`annotation_id`: uint8, `predicted_intensity`: float16),
/// addressed by the `(start, end)` ranges the pointer index above resolves.
/// Kept as one owned struct (not just two bare slices) so `Runner` can hold
/// it for the lifetime of a run and hand `&'db` slices to `Scorer`.
#[derive(Debug)]
pub struct FragmentIntensityArrays {
    _annotation_id_mmap: Mmap,
    _intensity_mmap: Mmap,
    pub annotation_id: &'static [u8],
    pub predicted_intensity: &'static [f16],
}

impl FragmentIntensityArrays {
    /// `dir` is the `arrays.mmappet` directory itself (a subdirectory of
    /// `git/featureprediction`'s `PredictionCache` -- this reader never
    /// touches that cache's `index.sqlite3`/`write.lock`, only its arrays).
    pub fn open(dir: &Path) -> Result<Self, FragmentIntensityCacheError> {
        let schema = read_mmappet_schema(dir)?;
        let i_intensity = mmappet_column_index(&schema, "predicted_intensity", "float16")?;
        let i_annotation = mmappet_column_index(&schema, "annotation_id", "uint8")?;

        // SAFETY: dtypes checked against schema.txt above.
        let (annotation_id_mmap, annotation_id) =
            unsafe { mmap_as_slice::<u8>(&dir.join(format!("{i_annotation}.bin")))? };
        let (intensity_mmap, predicted_intensity) =
            unsafe { mmap_as_slice::<f16>(&dir.join(format!("{i_intensity}.bin")))? };

        Ok(Self {
            _annotation_id_mmap: annotation_id_mmap,
            _intensity_mmap: intensity_mmap,
            annotation_id,
            predicted_intensity,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    /// Same manual-temp-dir convention as `pmsms.rs`'s tests -- no
    /// `tempfile` crate dependency in this workspace.
    fn temp_dir(name: &str) -> PathBuf {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "sage_fragment_intensity_cache_test_{}_{name}_{n}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create temp dir");
        dir
    }

    fn write_dataset(dir: &Path, annotation_id: &[u8], intensity: &[f16]) {
        std::fs::write(dir.join("schema.txt"), "float16 predicted_intensity\nuint8 annotation_id\n")
            .unwrap();
        let mut f0 = File::create(dir.join("0.bin")).unwrap();
        for v in intensity {
            f0.write_all(&v.to_le_bytes()).unwrap();
        }
        let mut f1 = File::create(dir.join("1.bin")).unwrap();
        f1.write_all(annotation_id).unwrap();
    }

    #[test]
    fn opens_and_slices_real_shaped_dataset() {
        let dir = temp_dir("opens_and_slices");
        let intensity = [f16::from_f32(0.1), f16::from_f32(0.5), f16::from_f32(1.0)];
        let annotation_id = [3u8, 90, 173];
        write_dataset(&dir, &annotation_id, &intensity);

        let arrays = FragmentIntensityArrays::open(&dir).unwrap();
        assert_eq!(arrays.annotation_id, &annotation_id);
        assert_eq!(arrays.predicted_intensity.len(), 3);
        assert!((arrays.predicted_intensity[1].to_f32() - 0.5).abs() < 1e-3);
    }

    #[test]
    fn rejects_wrong_dtype() {
        let dir = temp_dir("wrong_dtype");
        std::fs::write(dir.join("schema.txt"), "float32 predicted_intensity\nuint8 annotation_id\n")
            .unwrap();
        std::fs::write(dir.join("0.bin"), []).unwrap();
        std::fs::write(dir.join("1.bin"), []).unwrap();

        let err = FragmentIntensityArrays::open(&dir).unwrap_err();
        assert!(matches!(err, FragmentIntensityCacheError::WrongDtype(..)));
    }
}
