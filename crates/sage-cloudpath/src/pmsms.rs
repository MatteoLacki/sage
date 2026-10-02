//! Reader for the ionmaiden pmsms binary format.
//!
//! Inputs, given as explicit paths (never assumed to live together in one
//! directory with fixed filenames):
//!   pmsms.mmappet/        — schema.txt + one `<idx>.bin` file per column (see
//!                           git/mmappet's Python writer); columns looked up by
//!                           name. Needs `tof`(u32) and `intensity`(u32).
//!   tof2mz.mmappet/       — one float32 column `mz` indexed by tof; fragment m/z
//!                           is `mz[tof]`.
//!   <name>.mzcalib/       — optional fragment m/z recalibration (timstofu's
//!                           `MzRecalibration`): `bias.mmappet` plus an `mz/`
//!                           curve. When given, each fragment m/z is divided by
//!                           `1 + (bias + f_mz(mz) + fragment_shift_ppm)·1e-6`,
//!                           with `fragment_shift_ppm` read per precursor.
//!   precursors, either:
//!     .parquet            — read via the `parquet` crate
//!     .mmappet/           — schema.txt + one `<idx>.bin` file per column (see
//!                           git/mmappet's Python writer); columns looked up by name
//!   one row per precursor with columns: precursor_idx(u64), mz(f64), rt(f64),
//!   inv_ion_mobility(f64), charges(i64), fragment_spectrum_start(u64),
//!   fragment_event_cnt(u64)
//!   Two further columns, ppm_tol_lo(f64) and ppm_tol_hi(f64), are optional.
//!   When both are present in the schema, each row's `isolation_window` is
//!   set to `Tolerance::Ppm(ppm_tol_lo, ppm_tol_hi)`; when either is absent,
//!   `isolation_window` is `None` and the run's global `precursor_tol`
//!   applies instead — this keeps precursor tables without the columns
//!   working unchanged.
//!   A further column, intensity(i64), is also optional (see
//!   docs/ai/pmsms_input.md). So is fragment_shift_ppm(f64), required only
//!   with a recalibration.

#![cfg(feature = "parquet")]

use std::fs::File;
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use rayon::prelude::*;
use sage_core::mass::Tolerance;
use sage_core::spectrum::{
    Precursor, ProcessedSpectrum, RawSpectrum, Representation, SpectrumProcessor,
};

pub use crate::util::PmsmsPaths;

#[derive(Debug, thiserror::Error)]
pub enum PmsmsError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parquet error: {0}")]
    Parquet(#[from] parquet::errors::ParquetError),
    #[error("missing column: {0}")]
    MissingColumn(&'static str),
    #[error("mmappet column {0} has dtype `{1}`, expected `{2}`")]
    WrongDtype(&'static str, String, &'static str),
    #[error("tof {0} is outside the tof2mz table (length {1})")]
    TofOutOfRange(u32, usize),
    #[error("unrecognized precursors format (expected .parquet or .mmappet): {0}")]
    UnrecognizedPrecursorsFormat(PathBuf),
    #[error("fragment m/z recalibration {0} must have exactly one dimension, `mz`; found {1:?}")]
    MzRecalibrationDims(PathBuf, Vec<String>),
    #[error("fragment m/z recalibration {0}: {1}")]
    MzRecalibrationGrid(PathBuf, &'static str),
}

struct PrecursorRow {
    precursor_idx: u64,
    mz: f32,
    rt_minutes: f32,
    iim: f32,
    /// Decoded charge states (e.g. encoded value 234 → [2, 3, 4]).
    charges: Vec<u8>,
    frag_start: u64,
    frag_count: u64,
    /// Per-precursor ppm tolerance bounds, present only when both
    /// `ppm_tol_lo` and `ppm_tol_hi` columns exist in the precursors table.
    ppm_tol: Option<(f32, f32)>,
    ms1_intensity: Option<f32>,
    /// Per-precursor fragment m/z offset (ppm) added to the recalibration's
    /// `f_mz`; present only when the `fragment_shift_ppm` column exists.
    fragment_shift_ppm: Option<f64>,
}

/// Decode the pipeline's digit-concatenated charge encoding into individual charges.
/// e.g. 23 → [2, 3], 234 → [2, 3, 4], 2 → [2].
fn decode_charges(encoded: i64) -> Vec<u8> {
    if encoded <= 0 {
        return vec![];
    }
    let s = encoded.to_string();
    s.bytes()
        .filter_map(|b| {
            let d = b.wrapping_sub(b'0');
            if d > 0 && d <= 9 {
                Some(d)
            } else {
                None
            }
        })
        .collect()
}

fn read_precursors_parquet(parquet_path: &Path) -> Result<Vec<PrecursorRow>, PmsmsError> {
    let file = File::open(parquet_path)?;
    let reader = SerializedFileReader::new(file)?;

    // Locate column indices by name.
    let schema = reader.metadata().file_metadata().schema_descr();
    let col_idx = |name: &'static str| -> Result<usize, PmsmsError> {
        schema
            .columns()
            .iter()
            .position(|c| c.name() == name)
            .ok_or(PmsmsError::MissingColumn(name))
    };

    let i_mz = col_idx("mz")?;
    let i_rt = col_idx("rt")?;
    let i_iim = col_idx("inv_ion_mobility")?;
    let i_charge = col_idx("charges")?;
    let i_start = col_idx("fragment_spectrum_start")?;
    let i_count = col_idx("fragment_event_cnt")?;
    let i_pidx = col_idx("precursor_idx")?;
    // Optional: only present when the pipeline has computed a per-precursor
    // ppm tolerance. Absent columns mean every row falls back to the run's
    // global precursor_tol (see Precursor::effective_precursor_tol).
    let i_ppm_lo = col_idx("ppm_tol_lo").ok();
    let i_ppm_hi = col_idx("ppm_tol_hi").ok();
    let i_intensity = col_idx("intensity").ok();
    let i_fragment_shift = col_idx("fragment_shift_ppm").ok();

    let mut rows: Vec<PrecursorRow> = Vec::new();

    let iter = reader.get_row_iter(None)?;
    for row in iter {
        let row = row?;
        let precursor_idx = row.get_ulong(i_pidx)? as u64;
        let mz = row.get_double(i_mz)? as f32;
        let rt_minutes = (row.get_double(i_rt)? / 60.0) as f32;
        let iim = row.get_double(i_iim)? as f32;
        // charges column is int64 in the parquet; treat 0 as unknown
        let charge_raw = row.get_long(i_charge).unwrap_or(0);
        let charges = decode_charges(charge_raw);
        let frag_start = row.get_ulong(i_start)? as u64;
        let frag_count = row.get_ulong(i_count)? as u64;
        let ppm_tol = match (i_ppm_lo, i_ppm_hi) {
            (Some(lo), Some(hi)) => Some((row.get_double(lo)? as f32, row.get_double(hi)? as f32)),
            _ => None,
        };
        let ms1_intensity = match i_intensity {
            Some(i) => Some(row.get_long(i)? as f32),
            None => None,
        };
        let fragment_shift_ppm = match i_fragment_shift {
            Some(i) => Some(row.get_double(i)?),
            None => None,
        };

        if frag_count == 0 {
            continue;
        }

        rows.push(PrecursorRow {
            precursor_idx,
            mz,
            rt_minutes,
            iim,
            charges,
            frag_start,
            frag_count,
            ppm_tol,
            ms1_intensity,
            fragment_shift_ppm,
        });
    }

    Ok(rows)
}

/// Parse an mmappet `schema.txt` into an ordered list of (dtype, column name),
/// matching the numpy dtype strings written by git/mmappet's Python writer.
fn read_mmappet_schema(dir: &Path) -> Result<Vec<(String, String)>, PmsmsError> {
    let text = std::fs::read_to_string(dir.join("schema.txt"))?;
    Ok(text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| {
            let mut parts = line.splitn(2, ' ');
            let dtype = parts.next().unwrap_or_default().to_string();
            let name = parts.next().unwrap_or_default().trim().to_string();
            (dtype, name)
        })
        .collect())
}

/// Look up a named column's file index within an mmappet schema, asserting its
/// dtype matches what we're about to reinterpret its `<idx>.bin` file as.
fn mmappet_column_index(
    schema: &[(String, String)],
    name: &'static str,
    expected_dtype: &'static str,
) -> Result<usize, PmsmsError> {
    let (idx, (dtype, _)) = schema
        .iter()
        .enumerate()
        .find(|(_, (_, colname))| colname == name)
        .ok_or(PmsmsError::MissingColumn(name))?;
    if dtype != expected_dtype {
        return Err(PmsmsError::WrongDtype(name, dtype.clone(), expected_dtype));
    }
    Ok(idx)
}

/// Like [`mmappet_column_index`], but a missing column is `Ok(None)` rather
/// than an error — used for optional columns like `ppm_tol_lo`/`ppm_tol_hi`.
/// A present column with the wrong dtype still errors.
fn mmappet_optional_column_index(
    schema: &[(String, String)],
    name: &'static str,
    expected_dtype: &'static str,
) -> Result<Option<usize>, PmsmsError> {
    if !schema.iter().any(|(_, colname)| colname == name) {
        return Ok(None);
    }
    mmappet_column_index(schema, name, expected_dtype).map(Some)
}

fn read_precursors_mmappet(dir: &Path) -> Result<Vec<PrecursorRow>, PmsmsError> {
    let schema = read_mmappet_schema(dir)?;

    let i_pidx = mmappet_column_index(&schema, "precursor_idx", "uint64")?;
    let i_mz = mmappet_column_index(&schema, "mz", "float64")?;
    let i_rt = mmappet_column_index(&schema, "rt", "float64")?;
    let i_iim = mmappet_column_index(&schema, "inv_ion_mobility", "float64")?;
    let i_charge = mmappet_column_index(&schema, "charges", "int64")?;
    let i_start = mmappet_column_index(&schema, "fragment_spectrum_start", "uint64")?;
    let i_count = mmappet_column_index(&schema, "fragment_event_cnt", "uint64")?;
    let i_ppm_lo = mmappet_optional_column_index(&schema, "ppm_tol_lo", "float64")?;
    let i_ppm_hi = mmappet_optional_column_index(&schema, "ppm_tol_hi", "float64")?;
    let i_intensity = mmappet_optional_column_index(&schema, "intensity", "int64")?;
    let i_fragment_shift = mmappet_optional_column_index(&schema, "fragment_shift_ppm", "float64")?;

    // SAFETY: dtype checked against schema.txt above for each column.
    let (_m0, precursor_idx) = unsafe { mmap_as_slice::<u64>(&dir.join(format!("{i_pidx}.bin")))? };
    let (_m1, mz) = unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i_mz}.bin")))? };
    let (_m2, rt) = unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i_rt}.bin")))? };
    let (_m3, iim) = unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i_iim}.bin")))? };
    let (_m4, charges) = unsafe { mmap_as_slice::<i64>(&dir.join(format!("{i_charge}.bin")))? };
    let (_m5, frag_start) = unsafe { mmap_as_slice::<u64>(&dir.join(format!("{i_start}.bin")))? };
    let (_m6, frag_count) = unsafe { mmap_as_slice::<u64>(&dir.join(format!("{i_count}.bin")))? };
    // SAFETY: same as above — only mapped when the column exists in the schema.
    let ppm_lo = match i_ppm_lo {
        Some(i) => Some(unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i}.bin")))? }),
        None => None,
    };
    let ppm_hi = match i_ppm_hi {
        Some(i) => Some(unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i}.bin")))? }),
        None => None,
    };
    let intensity = match i_intensity {
        Some(i) => Some(unsafe { mmap_as_slice::<i64>(&dir.join(format!("{i}.bin")))? }),
        None => None,
    };
    let fragment_shift = match i_fragment_shift {
        Some(i) => Some(unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i}.bin")))? }),
        None => None,
    };

    let mut rows = Vec::with_capacity(precursor_idx.len());
    for i in 0..precursor_idx.len() {
        if frag_count[i] == 0 {
            continue;
        }
        let ppm_tol = match (&ppm_lo, &ppm_hi) {
            (Some((_, lo)), Some((_, hi))) => Some((lo[i] as f32, hi[i] as f32)),
            _ => None,
        };
        let ms1_intensity = intensity.as_ref().map(|(_, arr)| arr[i] as f32);
        let fragment_shift_ppm = fragment_shift.as_ref().map(|(_, arr)| arr[i]);
        rows.push(PrecursorRow {
            precursor_idx: precursor_idx[i],
            mz: mz[i] as f32,
            rt_minutes: (rt[i] / 60.0) as f32,
            iim: iim[i] as f32,
            charges: decode_charges(charges[i]),
            frag_start: frag_start[i],
            frag_count: frag_count[i],
            ppm_tol,
            ms1_intensity,
            fragment_shift_ppm,
        });
    }
    Ok(rows)
}

fn read_precursors(path: &Path) -> Result<Vec<PrecursorRow>, PmsmsError> {
    match path.extension().and_then(|ext| ext.to_str()) {
        Some("parquet") => read_precursors_parquet(path),
        Some("mmappet") => read_precursors_mmappet(path),
        _ => Err(PmsmsError::UnrecognizedPrecursorsFormat(path.to_path_buf())),
    }
}

/// Memory-map a binary file and reinterpret its contents as a typed slice.
///
/// # Safety
/// The file must contain valid data for type `T` (correct alignment, correct
/// element count). For our u32/f32 binary arrays this is always true when the
/// files were written by the ionmaiden pipeline.
unsafe fn mmap_as_slice<T: Copy>(path: &Path) -> Result<(Mmap, &'static [T]), PmsmsError> {
    let file = File::open(path)?;
    let mmap = Mmap::map(&file)?;
    let n = mmap.len() / std::mem::size_of::<T>();
    // SAFETY: the mmap outlives the slice through the returned owned Mmap.
    // The caller must hold the Mmap alive while the slice is used.
    let ptr = mmap.as_ptr() as *const T;
    let slice = std::slice::from_raw_parts(ptr, n);
    // Transmute to 'static lifetime — safe as long as Mmap is alive.
    // We return the Mmap with the slice so the caller owns both.
    let slice: &'static [T] = std::mem::transmute(slice);
    Ok((mmap, slice))
}

pub fn parse(paths: &PmsmsPaths, file_id: usize) -> Result<Vec<RawSpectrum>, PmsmsError> {
    parse_with(paths, file_id, |s| RawSpectrum {
        file_id: s.file_id,
        ms_level: s.ms_level,
        id: s.id,
        precursors: s.precursors,
        representation: s.representation,
        scan_start_time: s.scan_start_time,
        ion_injection_time: s.ion_injection_time,
        total_ion_current: s.total_ion_current,
        mz: s.mz.to_vec(),
        intensity: s.intensity.iter().map(|&i| i as f32).collect(),
        mobility: s.mobility,
    })
}

/// Process borrowed mapped fragment columns into owned, selected peak columns.
pub fn parse_processed(
    paths: &PmsmsPaths,
    file_id: usize,
    processor: &SpectrumProcessor,
) -> Result<Vec<ProcessedSpectrum>, PmsmsError> {
    parse_with(paths, file_id, |s| processor.process(s))
}

fn parse_with<T: Send>(
    paths: &PmsmsPaths,
    file_id: usize,
    transform: impl for<'a> Fn(RawSpectrum<&'a [f32], &'a [u32]>) -> T + Sync,
) -> Result<Vec<T>, PmsmsError> {
    let precursors = read_precursors(&paths.precursors)?;
    let recalibration = match &paths.mz_recalibration {
        Some(dir) => Some(FragmentMzRecalibration::load(dir)?),
        None => None,
    };

    let frag_schema = read_mmappet_schema(&paths.pmsms)?;
    if frag_schema.iter().any(|(_, name)| name == "mz") {
        log::warn!(
            "pmsms {} has an mz column: ignored, fragment m/z is read as tof2mz[tof]",
            paths.pmsms.display()
        );
    }
    let i_int = mmappet_column_index(&frag_schema, "intensity", "uint32")?;
    let i_tof = mmappet_column_index(&frag_schema, "tof", "uint32")?;
    let table_schema = read_mmappet_schema(&paths.tof2mz)?;
    let i_table = mmappet_column_index(&table_schema, "mz", "float32")?;
    // SAFETY: dtypes checked against schema.txt above.
    let (_frag_int_mmap, frag_int) =
        unsafe { mmap_as_slice::<u32>(&paths.pmsms.join(format!("{i_int}.bin")))? };
    let (_tof_mmap, tof) = unsafe { mmap_as_slice::<u32>(&paths.pmsms.join(format!("{i_tof}.bin")))? };
    let (_table_mmap, tof2mz) =
        unsafe { mmap_as_slice::<f32>(&paths.tof2mz.join(format!("{i_table}.bin")))? };

    let spectra = precursors
        .into_par_iter()
        .map(|p| {
            let start = p.frag_start as usize;
            let end = start + p.frag_count as usize;

            let mut mz = look_up_mz(&tof[start..end], tof2mz)?;
            if let Some(recalibration) = &recalibration {
                let shift_ppm = p
                    .fragment_shift_ppm
                    .ok_or(PmsmsError::MissingColumn("fragment_shift_ppm"))?;
                recalibration.correct(&mut mz, shift_ppm);
            }
            let int_slice = &frag_int[start..end];

            let total_ion_current = int_slice.iter().map(|&i| i as f32).sum();

            let isolation_window = p.ppm_tol.map(|(lo, hi)| Tolerance::Ppm(lo, hi));

            // One Precursor per charge state — mirrors how the MGF parser handles
            // CHARGE=234+ (regex extracts digits 2, 3, 4 as separate Precursor entries).
            let precursors: Vec<Precursor> = if p.charges.is_empty() {
                vec![Precursor {
                    mz: p.mz,
                    intensity: p.ms1_intensity,
                    charge: None,
                    spectrum_ref: None,
                    isolation_window,
                    inverse_ion_mobility: Some(p.iim),
                }]
            } else {
                p.charges
                    .iter()
                    .map(|&c| Precursor {
                        mz: p.mz,
                        intensity: p.ms1_intensity,
                        charge: Some(c),
                        spectrum_ref: None,
                        isolation_window,
                        inverse_ion_mobility: Some(p.iim),
                    })
                    .collect()
            };

            Ok(transform(RawSpectrum {
                file_id,
                ms_level: 2,
                id: format!("precursor_idx={}", p.precursor_idx),
                precursors,
                representation: Representation::Centroid,
                scan_start_time: p.rt_minutes,
                ion_injection_time: 0.0,
                total_ion_current,
                mz: &mz,
                intensity: int_slice,
                mobility: None,
            }))
        })
        .collect::<Result<Vec<T>, PmsmsError>>()?;

    Ok(spectra)
}

fn look_up_mz(tof: &[u32], tof2mz: &[f32]) -> Result<Vec<f32>, PmsmsError> {
    let mut mz = Vec::with_capacity(tof.len());
    for &t in tof {
        match tof2mz.get(t as usize) {
            Some(&value) => mz.push(value),
            None => return Err(PmsmsError::TofOutOfRange(t, tof2mz.len())),
        }
    }
    Ok(mz)
}

/// A piecewise-linear function on an evenly spaced grid, clamped to the edge
/// values outside it. Evaluates with the same f64 operations as timstofu's
/// `EvenlySpacedLinearSpline.njit_evaluator`, so SAGE and the Python export
/// path round recalibrated m/z to the same f32.
struct EvenlySpacedLinearSpline {
    x_min: f64,
    inv_spacing: f64,
    values: Vec<f64>,
}

impl EvenlySpacedLinearSpline {
    /// `values` must hold at least two points and `x_min < x_max`.
    fn new(x_min: f64, x_max: f64, values: Vec<f64>) -> Self {
        let spacing = (x_max - x_min) / (values.len() - 1) as f64;
        Self {
            x_min,
            inv_spacing: 1.0 / spacing,
            values,
        }
    }

    fn eval(&self, x: f64) -> f64 {
        let last = self.values.len() - 1;
        let idx_f = (x - self.x_min) * self.inv_spacing;
        if idx_f <= 0.0 {
            return self.values[0];
        }
        if idx_f >= last as f64 {
            return self.values[last];
        }
        let idx = idx_f as usize;
        let frac = idx_f - idx as f64;
        self.values[idx] * (1.0 - frac) + self.values[idx + 1] * frac
    }
}

/// A `.mzcalib` directory holding only the `mz` dimension: fragment ppm error
/// `bias + f_mz(mz)`, to which each precursor adds its own `fragment_shift_ppm`.
struct FragmentMzRecalibration {
    bias: f64,
    f_mz: EvenlySpacedLinearSpline,
}

impl FragmentMzRecalibration {
    fn load(dir: &Path) -> Result<Self, PmsmsError> {
        let dims = recalibration_dims(dir)?;
        if dims != ["mz"] {
            return Err(PmsmsError::MzRecalibrationDims(dir.to_path_buf(), dims));
        }
        let invalid = |reason| PmsmsError::MzRecalibrationGrid(dir.to_path_buf(), reason);

        let bias = read_f64_column(&dir.join("bias.mmappet"), "bias")?;
        let x_min = read_f64_column(&dir.join("mz/grid.mmappet"), "x_min")?;
        let x_max = read_f64_column(&dir.join("mz/grid.mmappet"), "x_max")?;
        let values = read_f64_column(&dir.join("mz/values.mmappet"), "y")?;
        if bias.len() != 1 || x_min.len() != 1 || x_max.len() != 1 {
            return Err(invalid("bias, x_min and x_max must each hold one value"));
        }
        if values.len() < 2 {
            return Err(invalid("the mz curve needs at least two grid points"));
        }
        if !(x_min[0] < x_max[0]) {
            return Err(invalid("the mz curve needs x_min < x_max"));
        }
        Ok(Self {
            bias: bias[0],
            f_mz: EvenlySpacedLinearSpline::new(x_min[0], x_max[0], values),
        })
    }

    /// `mz / (1 + (bias + f_mz(mz) + shift_ppm)·1e-6)` in f64, rounded to f32 —
    /// the arithmetic of timstofu's `apply_fragment_mz_recalibration_blocks`.
    fn correct(&self, mz: &mut [f32], shift_ppm: f64) {
        for value in mz.iter_mut() {
            let raw = *value as f64;
            let ppm = self.bias + self.f_mz.eval(raw) + shift_ppm;
            *value = (raw / (1.0 + ppm * 1e-6)) as f32;
        }
    }
}

/// Names of the curves in a `.mzcalib` directory: its subdirectories holding a
/// `grid.mmappet`, sorted (as timstofu's `MzRecalibration.load` finds them).
fn recalibration_dims(dir: &Path) -> Result<Vec<String>, PmsmsError> {
    let mut dims = Vec::new();
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.join("grid.mmappet").is_dir() {
            dims.push(path.file_name().unwrap_or_default().to_string_lossy().into_owned());
        }
    }
    dims.sort();
    Ok(dims)
}

fn read_f64_column(dir: &Path, name: &'static str) -> Result<Vec<f64>, PmsmsError> {
    let schema = read_mmappet_schema(dir)?;
    let i = mmappet_column_index(&schema, name, "float64")?;
    // SAFETY: dtype checked against schema.txt above.
    let (_map, values) = unsafe { mmap_as_slice::<f64>(&dir.join(format!("{i}.bin")))? };
    Ok(values.to_vec())
}

#[cfg(test)]
mod test {
    use super::*;

    fn fixture_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/pmsms_fixture")
    }

    fn pmsms_paths(pmsms: &Path, precursors: &Path, tof2mz: &Path, mz_recalibration: Option<&Path>) -> PmsmsPaths {
        PmsmsPaths {
            pmsms: pmsms.to_path_buf(),
            precursors: precursors.to_path_buf(),
            tof2mz: tof2mz.to_path_buf(),
            mz_recalibration: mz_recalibration.map(Path::to_path_buf),
        }
    }

    /// Dense tof2mz table rebuilt from the fixture's own `tof`/`mz` columns
    /// (its `mz` was once back-filled from a tof2mz table, so every tof maps to
    /// one m/z); tofs absent from the fixture map to 0.
    fn fixture_tof2mz_table() -> PathBuf {
        let pmsms_dir = fixture_dir().join("pmsms.mmappet");
        let schema = read_mmappet_schema(&pmsms_dir).unwrap();
        let i_tof = mmappet_column_index(&schema, "tof", "uint32").unwrap();
        let i_mz = mmappet_column_index(&schema, "mz", "float32").unwrap();
        let (_tof_map, tof) = unsafe { mmap_as_slice::<u32>(&pmsms_dir.join(format!("{i_tof}.bin"))).unwrap() };
        let (_mz_map, mz) = unsafe { mmap_as_slice::<f32>(&pmsms_dir.join(format!("{i_mz}.bin"))).unwrap() };

        let mut table = vec![0.0f32; *tof.iter().max().unwrap() as usize + 1];
        for (&t, &m) in tof.iter().zip(mz) {
            assert!(table[t as usize] == 0.0 || table[t as usize] == m, "tof {t} maps to two m/z");
            table[t as usize] = m;
        }
        MmappetDirBuilder::new("fixture_tof2mz").column("float32", "mz", &table).finish()
    }

    #[test]
    fn mmappet_and_parquet_precursors_agree() {
        let dir = fixture_dir();
        let pmsms_dir = dir.join("pmsms.mmappet");
        let table_dir = fixture_tof2mz_table();

        let via_parquet = parse(&pmsms_paths(&pmsms_dir, &dir.join("precursors.parquet"), &table_dir, None), 0)
            .expect("parse via parquet precursors");
        let via_mmappet = parse(&pmsms_paths(&pmsms_dir, &dir.join("precursors.mmappet"), &table_dir, None), 0)
            .expect("parse via mmappet precursors");

        assert!(!via_parquet.is_empty());
        assert_eq!(via_parquet.len(), via_mmappet.len());

        for (a, b) in via_parquet.iter().zip(via_mmappet.iter()) {
            assert_eq!(a.id, b.id);
            assert_eq!(a.mz, b.mz);
            assert_eq!(a.intensity, b.intensity);
            assert_eq!(a.precursors.len(), b.precursors.len());
            for (pa, pb) in a.precursors.iter().zip(b.precursors.iter()) {
                assert_eq!(pa.charge, pb.charge);
                assert!((pa.mz - pb.mz).abs() < 1e-6);
                assert_eq!(pa.inverse_ion_mobility, pb.inverse_ion_mobility);
                // Fixture has no ppm_tol_lo/ppm_tol_hi columns in either format —
                // both readers must fall back to `None` (run's global precursor_tol).
                assert!(pa.isolation_window.is_none());
                assert!(pb.isolation_window.is_none());
            }
        }
        let _ = std::fs::remove_dir_all(&table_dir);
    }

    #[test]
    fn mapped_processing_matches_owned_reader_for_both_precursor_formats() {
        let dir = fixture_dir();
        let table_dir = fixture_tof2mz_table();
        for name in ["precursors.parquet", "precursors.mmappet"] {
            let paths = pmsms_paths(&dir.join("pmsms.mmappet"), &dir.join(name), &table_dir, None);
            for deisotope in [false, true] {
                for top_n in [0, 2, usize::MAX] {
                    let processor = SpectrumProcessor::new(top_n, deisotope, 0.0, false);
                    let expected: Vec<_> = parse(&paths, 7)
                        .unwrap()
                        .into_iter()
                        .map(|s| processor.process(s))
                        .collect();
                    let actual = parse_processed(&paths, 7, &processor).unwrap();
                    // The reader's mappings are already closed here.
                    assert_eq!(actual.len(), expected.len());
                    for (a, b) in actual.iter().zip(&expected) {
                        assert_eq!(a.peaks, b.peaks);
                        assert_eq!(a.peak_charges, b.peak_charges);
                        assert_eq!(a.total_ion_current.to_bits(), b.total_ion_current.to_bits());
                        assert_eq!(a.id, b.id);
                        assert_eq!(a.file_id, b.file_id);
                        assert_eq!(a.level, b.level);
                        assert_eq!(a.scan_start_time, b.scan_start_time);
                        assert_eq!(a.ion_injection_time, b.ion_injection_time);
                        assert_eq!(a.precursors.len(), b.precursors.len());
                        for (pa, pb) in a.precursors.iter().zip(&b.precursors) {
                            assert_eq!(pa.mz, pb.mz);
                            assert_eq!(pa.charge, pb.charge);
                            assert_eq!(pa.intensity, pb.intensity);
                            assert_eq!(pa.spectrum_ref, pb.spectrum_ref);
                            assert_eq!(pa.inverse_ion_mobility, pb.inverse_ion_mobility);
                            assert_eq!(pa.isolation_window, pb.isolation_window);
                        }
                    }
                }
            }
        }
        let _ = std::fs::remove_dir_all(&table_dir);
    }

    #[test]
    fn unrecognized_precursors_format_errors() {
        let dir = fixture_dir();
        let table_dir = fixture_tof2mz_table();
        let bogus = dir.join("precursors.parquet").with_extension("txt");
        let err = parse(&pmsms_paths(&dir.join("pmsms.mmappet"), &bogus, &table_dir, None), 0).unwrap_err();
        assert!(matches!(err, PmsmsError::UnrecognizedPrecursorsFormat(_)));
        let _ = std::fs::remove_dir_all(&table_dir);
    }

    /// Build a minimal synthetic mmappet directory: `schema.txt` plus one
    /// `<idx>.bin` file per column, matching git/mmappet's on-disk format.
    struct MmappetDirBuilder {
        dir: PathBuf,
        schema_lines: Vec<String>,
        next_idx: usize,
    }

    impl MmappetDirBuilder {
        fn new(name: &str) -> Self {
            Self::at(temp_dir(name, "mmappet"))
        }

        /// A builder writing into `dir` (parents created), e.g. inside a `.mzcalib`.
        fn at(dir: PathBuf) -> Self {
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).expect("create mmappet dir");
            Self {
                dir,
                schema_lines: Vec::new(),
                next_idx: 0,
            }
        }

        /// Write one column: `dtype` must match git/mmappet's numpy dtype
        /// strings (e.g. "uint64", "float64", "int64", "float32", "uint32").
        fn column<T: Copy>(mut self, dtype: &str, name: &str, data: &[T]) -> Self {
            let idx = self.next_idx;
            self.next_idx += 1;
            self.schema_lines.push(format!("{dtype} {name}"));

            let mut bytes = Vec::with_capacity(data.len() * std::mem::size_of::<T>());
            for v in data {
                // SAFETY: T is Copy (plain numeric type); we just reinterpret
                // its bytes, matching how `mmap_as_slice` reads them back.
                let raw = unsafe {
                    std::slice::from_raw_parts(v as *const T as *const u8, std::mem::size_of::<T>())
                };
                bytes.extend_from_slice(raw);
            }
            std::fs::write(self.dir.join(format!("{idx}.bin")), bytes).expect("write column bin");
            self
        }

        fn finish(self) -> PathBuf {
            std::fs::write(self.dir.join("schema.txt"), self.schema_lines.join("\n"))
                .expect("write schema.txt");
            self.dir
        }
    }

    fn temp_dir(name: &str, extension: &str) -> PathBuf {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir().join(format!(
            "sage_pmsms_test_{}_{name}_{n}.{extension}",
            std::process::id()
        ))
    }

    /// 2 precursors x 2 fragments, tofs 1..=4, with a table mapping tof t to
    /// m/z 100·t; values are arbitrary filler for tests about precursor columns.
    fn mk_pmsms_fragments_dir(name: &str) -> (PathBuf, PathBuf) {
        let pmsms_dir = MmappetDirBuilder::new(name)
            .column("uint32", "tof", &[1u32, 2, 3, 4])
            .column("uint32", "intensity", &[10u32, 20, 30, 40])
            .finish();
        let table_dir = MmappetDirBuilder::new(&format!("{name}_tof2mz"))
            .column("float32", "mz", &[0.0f32, 100.0, 200.0, 300.0, 400.0])
            .finish();
        (pmsms_dir, table_dir)
    }

    #[test]
    fn ppm_tol_columns_present_set_isolation_window() {
        let precursors_dir = MmappetDirBuilder::new("precursors_with_ppm")
            .column("uint64", "precursor_idx", &[0u64, 1])
            .column("float64", "mz", &[500.0f64, 600.0])
            .column("float64", "rt", &[60.0f64, 120.0])
            .column("float64", "inv_ion_mobility", &[0.9f64, 1.0])
            .column("int64", "charges", &[2i64, 3])
            .column("uint64", "fragment_spectrum_start", &[0u64, 2])
            .column("uint64", "fragment_event_cnt", &[2u64, 2])
            .column("float64", "ppm_tol_lo", &[-10.0f64, -20.0])
            .column("float64", "ppm_tol_hi", &[10.0f64, 20.0])
            .finish();
        let (pmsms_dir, table_dir) = mk_pmsms_fragments_dir("pmsms_with_ppm");

        let spectra = parse(&pmsms_paths(&pmsms_dir, &precursors_dir, &table_dir, None), 0)
            .expect("parse with ppm columns");
        assert_eq!(spectra.len(), 2);

        match spectra[0].precursors[0].isolation_window {
            Some(Tolerance::Ppm(lo, hi)) => {
                assert_eq!(lo, -10.0);
                assert_eq!(hi, 10.0);
            }
            other => panic!("expected Some(Ppm(-10, 10)), got {other:?}"),
        }
        match spectra[1].precursors[0].isolation_window {
            Some(Tolerance::Ppm(lo, hi)) => {
                assert_eq!(lo, -20.0);
                assert_eq!(hi, 20.0);
            }
            other => panic!("expected Some(Ppm(-20, 20)), got {other:?}"),
        }

        for dir in [&precursors_dir, &pmsms_dir, &table_dir] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    #[test]
    fn ppm_tol_columns_absent_isolation_window_is_none() {
        // Same shape as above, but the ppm_tol_lo/ppm_tol_hi columns are
        // omitted entirely — this must parse exactly like today's real
        // pmsms tables (no error, isolation_window falls back to None).
        let precursors_dir = MmappetDirBuilder::new("precursors_without_ppm")
            .column("uint64", "precursor_idx", &[0u64, 1])
            .column("float64", "mz", &[500.0f64, 600.0])
            .column("float64", "rt", &[60.0f64, 120.0])
            .column("float64", "inv_ion_mobility", &[0.9f64, 1.0])
            .column("int64", "charges", &[2i64, 3])
            .column("uint64", "fragment_spectrum_start", &[0u64, 2])
            .column("uint64", "fragment_event_cnt", &[2u64, 2])
            .finish();
        let (pmsms_dir, table_dir) = mk_pmsms_fragments_dir("pmsms_without_ppm");

        let spectra = parse(&pmsms_paths(&pmsms_dir, &precursors_dir, &table_dir, None), 0)
            .expect("parse without ppm columns");
        assert_eq!(spectra.len(), 2);
        for spectrum in &spectra {
            assert!(spectrum.precursors[0].isolation_window.is_none());
        }

        for dir in [&precursors_dir, &pmsms_dir, &table_dir] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    fn two_precursors_builder(name: &str) -> MmappetDirBuilder {
        MmappetDirBuilder::new(name)
            .column("uint64", "precursor_idx", &[0u64, 1])
            .column("float64", "mz", &[500.0f64, 600.0])
            .column("float64", "rt", &[60.0f64, 120.0])
            .column("float64", "inv_ion_mobility", &[0.9f64, 1.0])
            .column("int64", "charges", &[2i64, 3])
            .column("uint64", "fragment_spectrum_start", &[0u64, 1])
            .column("uint64", "fragment_event_cnt", &[2u64, 2])
    }

    fn two_precursors_dir(name: &str) -> PathBuf {
        two_precursors_builder(name).finish()
    }

    #[test]
    fn tof2mz_table_supplies_fragment_mz_from_tof() {
        // Overlapping slices (0..2 and 1..3), as when precursors share a raw spectrum.
        let precursors_dir = two_precursors_dir("precursors_tof2mz");
        let events_dir = MmappetDirBuilder::new("events_tof_only")
            .column("uint32", "tof", &[1u32, 3, 4])
            .column("uint32", "intensity", &[10u32, 20, 30])
            .finish();
        let table_dir = MmappetDirBuilder::new("tof2mz_table")
            .column("float32", "mz", &[0.0f32, 101.5, 102.5, 103.5, 104.5])
            .finish();

        let spectra = parse(&pmsms_paths(&events_dir, &precursors_dir, &table_dir, None), 0)
            .expect("parse via tof2mz");
        assert_eq!(spectra[0].mz, vec![101.5f32, 103.5]);
        assert_eq!(spectra[1].mz, vec![103.5f32, 104.5]);
        assert_eq!(spectra[1].intensity, vec![20.0f32, 30.0]);

        for dir in [&precursors_dir, &events_dir, &table_dir] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    #[test]
    fn tof_outside_tof2mz_table_errors() {
        let precursors_dir = two_precursors_dir("precursors_tof_out_of_range");
        let events_dir = MmappetDirBuilder::new("events_tof_out_of_range")
            .column("uint32", "tof", &[1u32, 3, 9])
            .column("uint32", "intensity", &[10u32, 20, 30])
            .finish();
        let table_dir = MmappetDirBuilder::new("short_tof2mz_table")
            .column("float32", "mz", &[0.0f32, 101.5, 102.5, 103.5])
            .finish();

        let err = parse(&pmsms_paths(&events_dir, &precursors_dir, &table_dir, None), 0).unwrap_err();
        assert!(matches!(err, PmsmsError::TofOutOfRange(9, 4)));

        for dir in [&precursors_dir, &events_dir, &table_dir] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    // Generated with timstofu's numba code (necromerge2 venvs/common):
    // `EvenlySpacedLinearSpline(100.0, 1700.3, RECALIBRATION_CURVE).njit_evaluator()`
    // at each `RECALIBRATION_X`, and `apply_fragment_mz_recalibration_blocks`
    // over blocks [0, 4) and [4, 10) with shifts 1.5 and -2.25 ppm and bias 0.375.
    const RECALIBRATION_CURVE: [f64; 7] = [1.25, -0.5, 3.0, 2.2, -1.1, 0.4, 0.9];
    const RECALIBRATION_X_BITS: [u32; 10] = [
        0x42480000, 0x42c80000, 0x42c80001, 0x43b75bbc, 0x43d88fcc, 0x447a0000, 0x4491dccd, 0x44bb9000,
        0x44d4899a, 0x44e10000,
    ];
    const NUMBA_F_MZ_BITS: [u64; 10] = [
        0x3ff4000000000000,
        0x3ff4000000000000,
        0x3ff3fffff2900148,
        0xbfdfffff8d55603a,
        0x3fd7c56e02e95520,
        0x3feedde76c15d864,
        0xbff198d4e287c368,
        0x3fe0d072834f4db9,
        0x3feccccccccccccd,
        0x3feccccccccccccd,
    ];
    const NUMBA_CORRECTED_MZ_BITS: [u32; 10] = [
        0x4247ffd7, 0x42c7ffd7, 0x42c7ffd8, 0x43b75bab, 0x43d88fe1, 0x447a000f, 0x4491dce9, 0x44bb9011,
        0x44d489a8, 0x44e1000e,
    ];

    fn recalibration_x() -> Vec<f32> {
        RECALIBRATION_X_BITS.iter().map(|&b| f32::from_bits(b)).collect()
    }

    #[test]
    fn recalibration_matches_timstofu_bit_for_bit() {
        let recalibration = FragmentMzRecalibration {
            bias: 0.375,
            f_mz: EvenlySpacedLinearSpline::new(100.0, 1700.3, RECALIBRATION_CURVE.to_vec()),
        };
        let x = recalibration_x();
        for (&mz, &expected) in x.iter().zip(&NUMBA_F_MZ_BITS) {
            assert_eq!(recalibration.f_mz.eval(mz as f64).to_bits(), expected, "f_mz({mz})");
        }

        let mut corrected = x.clone();
        recalibration.correct(&mut corrected[..4], 1.5);
        recalibration.correct(&mut corrected[4..], -2.25);
        let corrected_bits: Vec<u32> = corrected.iter().map(|m| m.to_bits()).collect();
        assert_eq!(corrected_bits, NUMBA_CORRECTED_MZ_BITS);
    }

    /// A `.mzcalib` directory in timstofu's `MzRecalibration.dump` layout, with
    /// one curve per name in `dims`, each the `RECALIBRATION_CURVE` over [100, 1700.3].
    fn mk_mzcalib(name: &str, bias: f64, dims: &[&str]) -> PathBuf {
        let dir = temp_dir(name, "mzcalib");
        MmappetDirBuilder::at(dir.join("bias.mmappet"))
            .column("float64", "bias", &[bias])
            .finish();
        for dim in dims {
            MmappetDirBuilder::at(dir.join(dim).join("grid.mmappet"))
                .column("float64", "x_min", &[100.0f64])
                .column("float64", "x_max", &[1700.3f64])
                .finish();
            MmappetDirBuilder::at(dir.join(dim).join("values.mmappet"))
                .column("float64", "y", &RECALIBRATION_CURVE)
                .finish();
        }
        dir
    }

    /// Fragments are the ten `RECALIBRATION_X` m/z, at tofs 0..10; precursor 0
    /// owns rows 0..4 and precursor 1 rows 4..10.
    fn mk_recalibration_inputs(name: &str) -> (PathBuf, PathBuf) {
        let pmsms_dir = MmappetDirBuilder::new(name)
            .column("uint32", "tof", &(0u32..10).collect::<Vec<_>>())
            .column("uint32", "intensity", &[10u32; 10])
            .finish();
        let table_dir = MmappetDirBuilder::new(&format!("{name}_tof2mz"))
            .column("float32", "mz", &recalibration_x())
            .finish();
        (pmsms_dir, table_dir)
    }

    fn recalibration_precursors(name: &str) -> MmappetDirBuilder {
        MmappetDirBuilder::new(name)
            .column("uint64", "precursor_idx", &[0u64, 1])
            .column("float64", "mz", &[500.0f64, 600.0])
            .column("float64", "rt", &[60.0f64, 120.0])
            .column("float64", "inv_ion_mobility", &[0.9f64, 1.0])
            .column("int64", "charges", &[2i64, 3])
            .column("uint64", "fragment_spectrum_start", &[0u64, 4])
            .column("uint64", "fragment_event_cnt", &[4u64, 6])
    }

    #[test]
    fn mz_recalibration_applies_curve_bias_and_each_precursors_shift() {
        let (pmsms_dir, table_dir) = mk_recalibration_inputs("pmsms_recalibrated");
        let precursors_dir = recalibration_precursors("precursors_with_shift")
            .column("float64", "fragment_shift_ppm", &[1.5f64, -2.25])
            .finish();
        let mzcalib = mk_mzcalib("mz_only", 0.375, &["mz"]);

        let spectra = parse(&pmsms_paths(&pmsms_dir, &precursors_dir, &table_dir, Some(&mzcalib)), 0)
            .expect("parse with recalibration");
        let corrected_bits: Vec<u32> = spectra.iter().flat_map(|s| s.mz.iter().map(|m| m.to_bits())).collect();
        assert_eq!(corrected_bits, NUMBA_CORRECTED_MZ_BITS);

        for dir in [&pmsms_dir, &table_dir, &precursors_dir, &mzcalib] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    #[test]
    fn mz_recalibration_without_fragment_shift_column_errors() {
        let (pmsms_dir, table_dir) = mk_recalibration_inputs("pmsms_no_shift");
        let precursors_dir = recalibration_precursors("precursors_without_shift").finish();
        let mzcalib = mk_mzcalib("mz_only_no_shift", 0.0, &["mz"]);

        let err = parse(&pmsms_paths(&pmsms_dir, &precursors_dir, &table_dir, Some(&mzcalib)), 0).unwrap_err();
        assert!(matches!(err, PmsmsError::MissingColumn("fragment_shift_ppm")));

        for dir in [&pmsms_dir, &table_dir, &precursors_dir, &mzcalib] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }

    #[test]
    fn mz_recalibration_needs_exactly_the_mz_dimension() {
        let (pmsms_dir, table_dir) = mk_recalibration_inputs("pmsms_bad_dims");
        let precursors_dir = recalibration_precursors("precursors_bad_dims")
            .column("float64", "fragment_shift_ppm", &[0.0f64, 0.0])
            .finish();
        for dims in [&["mz", "rt"][..], &["rt"][..]] {
            let mzcalib = mk_mzcalib("bad_dims", 0.0, dims);
            let err =
                parse(&pmsms_paths(&pmsms_dir, &precursors_dir, &table_dir, Some(&mzcalib)), 0).unwrap_err();
            match err {
                PmsmsError::MzRecalibrationDims(_, found) => assert_eq!(found, dims),
                other => panic!("expected MzRecalibrationDims, got {other:?}"),
            }
            let _ = std::fs::remove_dir_all(&mzcalib);
        }

        for dir in [&pmsms_dir, &table_dir, &precursors_dir] {
            let _ = std::fs::remove_dir_all(dir);
        }
    }
}
