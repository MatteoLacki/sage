//! Reader for a SPECLIB's predicted fragment intensities (format version 3), as
//! written by `git/featureprediction`'s `feature-prediction-make-speclib` or
//! `predict_prosit` (see `docs/ai/predicted_fragment_intensity.md`).
//!
//! Two separate inputs, given as explicit paths:
//!   - the index (`index.mmappet`), positional: one `start` column (uint64,
//!     one per (peptide row, charge) slot plus the end), slot `i` being rows
//!     `start[i]..start[i + 1]` of the arrays below, an empty range meaning no
//!     prediction; slots numbered by the SPECLIB's `slot_order`
//!     (`predicted_properties::read_slot_order`). Its `meta.json` carries a
//!     `dumped_peptides_sha256` the caller must check (see
//!     `dumped_peptides_fingerprint`).
//!   - the `arrays.mmappet/` directory (schema.txt +
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
use sage_core::intensity_codec::FragmentIntensities;

use crate::predicted_properties::{
    read_fingerprint, read_speclib_manifest, MappedColumn, PredictedPropertiesError,
};

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
    #[error("invalid fragment intensity format: {0}")]
    InvalidFormat(String),
    #[error("invalid fragment intensity metadata: {0}")]
    Json(#[from] serde_json::Error),
}

/// A SPECLIB's index, mapped -- see module docs.
#[derive(Clone)]
pub struct FragmentIntensityIndex {
    pub start: MappedColumn<u64>,
}

/// The `meta.json` equations this reader implements.
const FRAGMENTS_AT_POSITION: &str = "arrays[start[position] : start[position + 1]]";
const SQRT8: &str = "round_half_even(255 * sqrt(intensity / max(intensity at position)))";
const FP16: &str = "intensity as float16";

/// Map the index in `dir`. Returns the `dumped_peptides` fingerprint alongside it
/// -- the caller must check it against the current run's own digest before trusting
/// these positions.
pub fn read_fragment_intensity_index(
    dir: &Path,
) -> Result<(String, FragmentIntensityIndex), FragmentIntensityCacheError> {
    let meta = read_speclib_manifest(dir)?;
    let equation = meta.get("fragments_at_position").and_then(|v| v.as_str());
    if equation != Some(FRAGMENTS_AT_POSITION) {
        return Err(FragmentIntensityCacheError::InvalidFormat(format!(
            "{} has fragments_at_position {equation:?}, which this SAGE does not read",
            dir.display()
        )));
    }
    let fingerprint = read_fingerprint(dir)?;
    let start = MappedColumn::open(dir, "start")?;
    Ok((fingerprint, FragmentIntensityIndex { start }))
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
    if file.metadata()?.len() % std::mem::size_of::<T>() as u64 != 0 {
        return Err(FragmentIntensityCacheError::InvalidFormat(format!("{} has a partial element", path.display())));
    }
    let mmap = Mmap::map(&file)?;
    let n = mmap.len() / std::mem::size_of::<T>();
    let ptr = mmap.as_ptr() as *const T;
    let slice = std::slice::from_raw_parts(ptr, n);
    Ok((mmap, slice))
}

/// The shared, ever-growing `arrays.mmappet` dataset -- two flat mmapped
/// columns (`annotation_id`: uint8, `predicted_intensity`: float16 or uint8, as
/// its `meta.json`'s `intensity_encoding` declares: `fp16` or `sqrt8_basepeak_v1`),
/// addressed by the `(start, end)` ranges the pointer index above resolves.
/// Kept as one owned struct (not just two bare slices) so `Runner` can hold
/// it for the lifetime of a run and hand `&'db` slices to `Scorer`.
#[derive(Debug)]
pub struct FragmentIntensityArrays {
    _annotation_id_mmap: Mmap,
    _intensity_mmap: Mmap,
    pub annotation_id: &'static [u8],
    pub predicted_intensity: FragmentIntensities<'static>,
}

impl FragmentIntensityArrays {
    /// `dir` is the SPECLIB's `arrays.mmappet` directory itself.
    pub fn open(dir: &Path) -> Result<Self, FragmentIntensityCacheError> {
        let schema = read_mmappet_schema(dir)?;
        let manifest = read_speclib_manifest(dir)?;
        let sqrt8 = match manifest.get("predicted_intensity").and_then(|v| v.as_str()) {
            Some(SQRT8) => true,
            Some(FP16) => false,
            other => {
                return Err(FragmentIntensityCacheError::InvalidFormat(format!(
                    "{} has predicted_intensity {other:?}, which this SAGE does not read",
                    dir.display()
                )))
            }
        };
        let dtype = if sqrt8 { "uint8" } else { "float16" };
        let i_intensity = mmappet_column_index(&schema, "predicted_intensity", dtype)?;
        let i_annotation = mmappet_column_index(&schema, "annotation_id", "uint8")?;

        // SAFETY: dtypes checked against schema.txt above.
        let (annotation_id_mmap, annotation_id) =
            unsafe { mmap_as_slice::<u8>(&dir.join(format!("{i_annotation}.bin")))? };
        let (intensity_mmap, predicted_intensity) = if sqrt8 {
            let (map, codes) = unsafe { mmap_as_slice::<u8>(&dir.join(format!("{i_intensity}.bin")))? };
            (map, FragmentIntensities::Sqrt8Basepeak(codes))
        } else {
            let (map, values) = unsafe { mmap_as_slice::<f16>(&dir.join(format!("{i_intensity}.bin")))? };
            (map, FragmentIntensities::Fp16(values))
        };
        if annotation_id.len() != predicted_intensity.len() {
            return Err(FragmentIntensityCacheError::InvalidFormat("intensities and annotation IDs differ in length".into()));
        }

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

    /// A format-4 fragments part's meta.json with `predicted_intensity` = `equation`.
    fn write_manifest(dir: &Path, equation: &str) {
        let meta = serde_json::json!({
            "format": "speclib", "format_version": 4, "part": "fragments",
            "predicted_intensity": equation,
        });
        std::fs::write(dir.join("meta.json"), meta.to_string()).unwrap();
    }

    fn write_dataset(dir: &Path, annotation_id: &[u8], intensity: &[f16]) {
        std::fs::write(dir.join("schema.txt"), "float16 predicted_intensity\nuint8 annotation_id\n")
            .unwrap();
        write_manifest(dir, FP16);
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
        assert!((arrays.predicted_intensity.value(1) - 0.5).abs() < 1e-3);
    }

    #[test]
    fn rejects_wrong_dtype() {
        let dir = temp_dir("wrong_dtype");
        std::fs::write(dir.join("schema.txt"), "float32 predicted_intensity\nuint8 annotation_id\n")
            .unwrap();
        write_manifest(&dir, FP16);
        std::fs::write(dir.join("0.bin"), []).unwrap();
        std::fs::write(dir.join("1.bin"), []).unwrap();

        let err = FragmentIntensityArrays::open(&dir).unwrap_err();
        assert!(matches!(err, FragmentIntensityCacheError::WrongDtype(..)));
    }

    #[test]
    fn requires_a_declared_intensity_encoding() {
        let dir = temp_dir("undeclared");
        write_dataset(&dir, &[3u8], &[f16::from_f32(1.0)]);
        std::fs::write(dir.join("meta.json"), r#"{"format":"speclib","format_version":3}"#).unwrap();
        let err = FragmentIntensityArrays::open(&dir).unwrap_err();
        assert!(err.to_string().contains("not part of a SPECLIB of format version 4"), "{err}");

        write_manifest(&dir, "fp16");
        let err = FragmentIntensityArrays::open(&dir).unwrap_err();
        assert!(err.to_string().contains("predicted_intensity Some(\"fp16\")"), "{err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn sqrt8_requires_explicit_metadata_and_matching_columns() {
        let dir = temp_dir("sqrt8");
        std::fs::write(dir.join("schema.txt"), "uint8 predicted_intensity\nuint8 annotation_id\n").unwrap();
        std::fs::write(dir.join("0.bin"), [0, 128, 255]).unwrap();
        std::fs::write(dir.join("1.bin"), [0, 90, 173]).unwrap();
        assert!(FragmentIntensityArrays::open(&dir).is_err());
        write_manifest(&dir, "sqrt8_basepeak_v1");
        assert!(FragmentIntensityArrays::open(&dir).is_err());
        write_manifest(&dir, SQRT8);
        {
            let arrays = FragmentIntensityArrays::open(&dir).unwrap();
            assert_eq!(arrays.predicted_intensity.codes().unwrap(), [0, 128, 255]);
            assert_eq!(arrays.predicted_intensity.value(0), 0.0);
            assert_eq!(arrays.predicted_intensity.value(2), 1.0);
        }
        std::fs::write(dir.join("1.bin"), [0, 90]).unwrap();
        assert!(FragmentIntensityArrays::open(&dir).is_err());
        std::fs::remove_dir_all(dir).unwrap();
    }
}
