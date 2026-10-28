//! Readers for externally-predicted peptide RT/IIM properties (see
//! `necromerge2`'s `git/featureprediction` and `plans/speclib_rt_iim.md`). RT and
//! IIM are two independent, separately-loadable files -- RT has no charge dimension
//! (Chronologer_RT doesn't take one), IIM does.
//!
//! Both are positional, not `sequence`-keyed: row `i` of the per-run calibrated
//! `predicted_rt.mmappet` (float64 minutes) is peptide row `i` of the
//! `dumped_peptides` it was derived from; IIM is a SPECLIB's `iim.mmappet`
//! (float32 1/K0), one value per (peptide row, charge) slot in the SPECLIB's
//! `slot_order` (`read_slot_order`). Each is an mmappet dataset (`schema.txt` +
//! `0.bin`) mapped in place, with a `meta.json` carrying `dumped_peptides_sha256`
//! (see `dumped_peptides_fingerprint`) -- the caller (`sage-cli`'s `Runner`) must
//! check it against the current run's own digest before trusting these positions.
//! A per-charge file's `meta.json` also carries the `min_charge`/`max_charge` its
//! slots cover (`read_charge_range`).

use memmap2::Mmap;
use sage_core::predicted::SlotOrder;
use std::marker::PhantomData;
use std::path::Path;
use std::sync::Arc;

#[derive(Debug, thiserror::Error)]
pub enum PredictedPropertiesError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("{0} has no column `{1}`")]
    MissingColumn(String, &'static str),
    #[error("{0}: column `{1}` has dtype `{2}`, expected `{3}`")]
    WrongDtype(String, &'static str, String, &'static str),
    #[error("{0}: column `{1}` is not a whole number of `{2}` values")]
    BadLength(String, &'static str, &'static str),
    #[error("{0} is missing `dumped_peptides_sha256` in its meta.json -- regenerate it with a current `feature-prediction` build")]
    MissingFingerprint(String),
    #[error("{0} is missing `min_charge`/`max_charge` in its meta.json -- regenerate it with a current `feature-prediction` or `predict_prosit` build")]
    MissingChargeRange(String),
    #[error("{0} is not part of a SPECLIB of format version {SPECLIB_FORMAT_VERSION} (its meta.json says format {1}, version {2}) -- regenerate it with `feature-prediction-make-speclib` or `predict_prosit`")]
    NotSpeclib(String, String, String),
    #[error("{0} has position_of_row_and_charge {1}, which this SAGE does not read")]
    UnknownSlotOrder(String, String),
}

/// The SPECLIB format version this SAGE reads.
pub const SPECLIB_FORMAT_VERSION: u64 = 4;

/// Element types a mapped column can hold, by their numpy dtype name.
pub trait ColumnType: Copy + Send + Sync + 'static {
    const DTYPE: &'static str;
}

impl ColumnType for f64 {
    const DTYPE: &'static str = "float64";
}

impl ColumnType for i64 {
    const DTYPE: &'static str = "int64";
}

impl ColumnType for f32 {
    const DTYPE: &'static str = "float32";
}

impl ColumnType for u64 {
    const DTYPE: &'static str = "uint64";
}

/// One column of an mmappet dataset, mapped read-only; clones share the mapping.
#[derive(Clone)]
pub struct MappedColumn<T: ColumnType> {
    map: Option<Arc<Mmap>>,
    len: usize,
    element: PhantomData<T>,
}

impl<T: ColumnType> MappedColumn<T> {
    /// Column `name` of the mmappet dataset `dir`, which must hold `T`.
    pub fn open(dir: &Path, name: &'static str) -> Result<Self, PredictedPropertiesError> {
        let label = dir.display().to_string();
        let schema = std::fs::read_to_string(dir.join("schema.txt"))?;
        let (index, dtype) = schema
            .lines()
            .filter_map(|line| line.split_once(' '))
            .enumerate()
            .find(|(_, (_, column))| column.trim() == name)
            .map(|(index, (dtype, _))| (index, dtype.to_string()))
            .ok_or_else(|| PredictedPropertiesError::MissingColumn(label.clone(), name))?;
        if dtype != T::DTYPE {
            return Err(PredictedPropertiesError::WrongDtype(
                label,
                name,
                dtype,
                T::DTYPE,
            ));
        }
        let file = std::fs::File::open(dir.join(format!("{index}.bin")))?;
        let bytes = file.metadata()?.len() as usize;
        if bytes % std::mem::size_of::<T>() != 0 {
            return Err(PredictedPropertiesError::BadLength(label, name, T::DTYPE));
        }
        let map = match bytes {
            0 => None,
            // SAFETY: prediction files are written once and never modified while mapped.
            _ => Some(Arc::new(unsafe { Mmap::map(&file)? })),
        };
        Ok(Self {
            map,
            len: bytes / std::mem::size_of::<T>(),
            element: PhantomData,
        })
    }

    pub fn as_slice(&self) -> &[T] {
        match &self.map {
            None => &[],
            // SAFETY: the mapping starts page-aligned, holds `len` values of `T` (checked
            // in `open`), every bit pattern is a valid `T`, and it lives as long as `self`.
            Some(map) => unsafe { std::slice::from_raw_parts(map.as_ptr() as *const T, self.len) },
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// The `dumped_peptides_sha256` in `dir/meta.json`.
pub fn read_fingerprint(dir: &Path) -> Result<String, PredictedPropertiesError> {
    let meta: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("meta.json"))?)?;
    meta.get("dumped_peptides_sha256")
        .and_then(|value| value.as_str())
        .map(str::to_string)
        .ok_or_else(|| PredictedPropertiesError::MissingFingerprint(dir.display().to_string()))
}

/// The `(min_charge, max_charge)` in `dir/meta.json`: the precursor charges the slots
/// of a per-charge file (dense over `(peptide_row, charge)`) cover.
pub fn read_charge_range(dir: &Path) -> Result<(u8, u8), PredictedPropertiesError> {
    let meta: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("meta.json"))?)?;
    match (meta_charge(&meta, "min_charge"), meta_charge(&meta, "max_charge")) {
        (Some(min_charge), Some(max_charge)) => Ok((min_charge, max_charge)),
        _ => Err(PredictedPropertiesError::MissingChargeRange(dir.display().to_string())),
    }
}

fn meta_charge(meta: &serde_json::Value, key: &str) -> Option<u8> {
    let charge = meta.get(key)?.as_u64()?;
    u8::try_from(charge).ok()
}

/// Column `rt` (float64, minutes), one row per `dumped_peptides` row, in that same
/// order. Returns the `dumped_peptides` fingerprint alongside the mapped values.
pub fn read_predicted_rt(
    dir: &Path,
) -> Result<(String, MappedColumn<f64>), PredictedPropertiesError> {
    Ok((read_fingerprint(dir)?, MappedColumn::open(dir, "rt")?))
}

/// A SPECLIB dataset's `meta.json` (its manifest), refused unless it declares
/// SPECLIB format version [`SPECLIB_FORMAT_VERSION`].
pub fn read_speclib_manifest(dir: &Path) -> Result<serde_json::Value, PredictedPropertiesError> {
    let meta: serde_json::Value = serde_json::from_slice(&std::fs::read(dir.join("meta.json"))?)?;
    let format = meta.get("format").and_then(|value| value.as_str());
    let version = meta.get("format_version").and_then(|value| value.as_u64());
    if format != Some("speclib") || version != Some(SPECLIB_FORMAT_VERSION) {
        return Err(PredictedPropertiesError::NotSpeclib(
            dir.display().to_string(),
            format!("{format:?}"),
            format!("{version:?}"),
        ));
    }
    Ok(meta)
}

/// Where the SPECLIB dataset in `dir` keeps each (peptide row, charge).
pub fn read_slot_order(dir: &Path) -> Result<SlotOrder, PredictedPropertiesError> {
    let meta = read_speclib_manifest(dir)?;
    let equation = meta.get("position_of_row_and_charge").and_then(|value| value.as_str());
    equation.and_then(SlotOrder::from_equation).ok_or_else(|| {
        PredictedPropertiesError::UnknownSlotOrder(dir.display().to_string(), format!("{equation:?}"))
    })
}

/// Column `iim` (float32, 1/K0) of a SPECLIB's `iim.mmappet`, one value per slot
/// (`read_slot_order`). Returns the `dumped_peptides` fingerprint alongside the
/// mapped values.
pub fn read_predicted_iim(
    dir: &Path,
) -> Result<(String, MappedColumn<f32>), PredictedPropertiesError> {
    read_speclib_manifest(dir)?;
    Ok((read_fingerprint(dir)?, MappedColumn::open(dir, "iim")?))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn temp_dir(name: &str) -> PathBuf {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "sage_predicted_properties_test_{}_{name}_{n}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    fn write_rt(dir: &Path, dtype: &str, values: &[f64]) {
        std::fs::write(dir.join("schema.txt"), format!("{dtype} rt")).unwrap();
        let bytes: Vec<u8> = values.iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(dir.join("0.bin"), bytes).unwrap();
        std::fs::write(
            dir.join("meta.json"),
            r#"{"dumped_peptides_sha256": "abc"}"#,
        )
        .unwrap();
    }

    #[test]
    fn maps_values_and_fingerprint() {
        let dir = temp_dir("maps");
        write_rt(&dir, "float64", &[1.5, f64::NAN, 30.25]);
        let (fingerprint, rt) = read_predicted_rt(&dir).unwrap();
        assert_eq!(fingerprint, "abc");
        assert_eq!(rt.len(), 3);
        assert_eq!(rt.as_slice()[0], 1.5);
        assert!(rt.as_slice()[1].is_nan());
        assert_eq!(rt.clone().as_slice()[2], 30.25);
    }

    #[test]
    fn rejects_wrong_dtype_missing_column_and_fingerprint() {
        let dir = temp_dir("rejects");
        write_rt(&dir, "float32", &[1.0]);
        assert!(matches!(
            read_predicted_rt(&dir),
            Err(PredictedPropertiesError::WrongDtype(..))
        ));
        write_rt(&dir, "float64", &[1.0]);
        std::fs::write(dir.join("meta.json"), "{}").unwrap();
        assert!(matches!(
            read_predicted_rt(&dir),
            Err(PredictedPropertiesError::MissingFingerprint(..))
        ));
    }

    #[test]
    fn reads_charge_range_and_rejects_its_absence() {
        let dir = temp_dir("charges");
        write_rt(&dir, "float64", &[1.0]);
        assert!(matches!(
            read_charge_range(&dir),
            Err(PredictedPropertiesError::MissingChargeRange(..))
        ));
        std::fs::write(
            dir.join("meta.json"),
            r#"{"dumped_peptides_sha256": "abc", "min_charge": 2, "max_charge": 4}"#,
        )
        .unwrap();
        assert_eq!(read_charge_range(&dir).unwrap(), (2, 4));
    }

    fn write_iim(dir: &Path, dtype: &str, meta: &str) {
        std::fs::write(dir.join("schema.txt"), format!("{dtype} iim")).unwrap();
        let bytes: Vec<u8> = [1.25f32, f32::NAN].iter().flat_map(|v| v.to_le_bytes()).collect();
        std::fs::write(dir.join("0.bin"), bytes).unwrap();
        std::fs::write(dir.join("meta.json"), meta).unwrap();
    }

    const IIM_META: &str = r#"{"format": "speclib", "format_version": 4, "position_of_row_and_charge": "(charge - min_charge) * row_count + row", "dumped_peptides_sha256": "abc", "min_charge": 2, "max_charge": 2}"#;

    #[test]
    fn maps_speclib_iim_and_its_slot_order() {
        let dir = temp_dir("iim");
        write_iim(&dir, "float32", IIM_META);
        let (fingerprint, iim) = read_predicted_iim(&dir).unwrap();
        assert_eq!(fingerprint, "abc");
        assert_eq!(iim.as_slice()[0], 1.25);
        assert!(iim.as_slice()[1].is_nan());
        assert_eq!(read_slot_order(&dir).unwrap(), SlotOrder::ChargeMajor);
    }

    #[test]
    fn refuses_iim_outside_speclib_format_4() {
        let dir = temp_dir("iim_refused");
        write_iim(&dir, "float64", IIM_META);
        assert!(matches!(
            read_predicted_iim(&dir),
            Err(PredictedPropertiesError::WrongDtype(..))
        ));
        write_iim(&dir, "float32", &IIM_META.replace(r#""format_version": 4"#, r#""format_version": 3"#));
        assert!(matches!(
            read_predicted_iim(&dir),
            Err(PredictedPropertiesError::NotSpeclib(..))
        ));
        write_iim(&dir, "float32", &IIM_META.replace("(charge - min_charge) * row_count + row", "row"));
        assert!(matches!(
            read_slot_order(&dir),
            Err(PredictedPropertiesError::UnknownSlotOrder(..))
        ));
    }
}
