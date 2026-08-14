//! Reader for the ionmaiden pmsms binary format.
//!
//! Two inputs, given as explicit paths (never assumed to live together in one
//! directory with fixed filenames):
//!   pmsms.mmappet/        — schema.txt + one `<idx>.bin` file per column (see
//!                           git/mmappet's Python writer); columns looked up by
//!                           name. Needs `mz`(f32) and `intensity`(u32).
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

#![cfg(feature = "parquet")]

use std::fs::File;
use std::path::{Path, PathBuf};

use memmap2::Mmap;
use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;
use sage_core::mass::Tolerance;
use sage_core::spectrum::{Precursor, RawSpectrum, Representation};

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
    #[error("unrecognized precursors format (expected .parquet or .mmappet): {0}")]
    UnrecognizedPrecursorsFormat(PathBuf),
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
            if d > 0 && d <= 9 { Some(d) } else { None }
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
        return Err(PmsmsError::WrongDtype(
            name,
            dtype.clone(),
            expected_dtype,
        ));
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

    let mut rows = Vec::with_capacity(precursor_idx.len());
    for i in 0..precursor_idx.len() {
        if frag_count[i] == 0 {
            continue;
        }
        let ppm_tol = match (&ppm_lo, &ppm_hi) {
            (Some((_, lo)), Some((_, hi))) => Some((lo[i] as f32, hi[i] as f32)),
            _ => None,
        };
        rows.push(PrecursorRow {
            precursor_idx: precursor_idx[i],
            mz: mz[i] as f32,
            rt_minutes: (rt[i] / 60.0) as f32,
            iim: iim[i] as f32,
            charges: decode_charges(charges[i]),
            frag_start: frag_start[i],
            frag_count: frag_count[i],
            ppm_tol,
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

pub fn parse(
    pmsms_dir: &Path,
    precursors_path: &Path,
    file_id: usize,
) -> Result<Vec<RawSpectrum>, PmsmsError> {
    let precursors = read_precursors(precursors_path)?;

    let frag_schema = read_mmappet_schema(pmsms_dir)?;
    let i_mz = mmappet_column_index(&frag_schema, "mz", "float32")?;
    let i_int = mmappet_column_index(&frag_schema, "intensity", "uint32")?;

    // Memory-map both fragment column arrays.
    // SAFETY: dtype checked against schema.txt above for each column.
    let (_frag_mz_mmap, frag_mz) =
        unsafe { mmap_as_slice::<f32>(&pmsms_dir.join(format!("{i_mz}.bin")))? };
    let (_frag_int_mmap, frag_int) =
        unsafe { mmap_as_slice::<u32>(&pmsms_dir.join(format!("{i_int}.bin")))? };

    let mut spectra = Vec::with_capacity(precursors.len());

    for p in &precursors {
        let start = p.frag_start as usize;
        let end = start + p.frag_count as usize;

        let mz_vec = frag_mz[start..end].to_vec();
        let int_slice = &frag_int[start..end];

        let int_vec: Vec<f32> = int_slice.iter().map(|&i| i as f32).collect();
        let total_ion_current: f32 = int_vec.iter().sum();

        let isolation_window = p.ppm_tol.map(|(lo, hi)| Tolerance::Ppm(lo, hi));

        // One Precursor per charge state — mirrors how the MGF parser handles
        // CHARGE=234+ (regex extracts digits 2, 3, 4 as separate Precursor entries).
        let precursors: Vec<Precursor> = if p.charges.is_empty() {
            vec![Precursor {
                mz: p.mz,
                intensity: None,
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
                    intensity: None,
                    charge: Some(c),
                    spectrum_ref: None,
                    isolation_window,
                    inverse_ion_mobility: Some(p.iim),
                })
                .collect()
        };

        spectra.push(RawSpectrum {
            file_id,
            ms_level: 2,
            id: format!("precursor_idx={}", p.precursor_idx),
            precursors,
            representation: Representation::Centroid,
            scan_start_time: p.rt_minutes,
            ion_injection_time: 0.0,
            total_ion_current,
            mz: mz_vec,
            intensity: int_vec,
            mobility: None,
        });
    }

    Ok(spectra)
}

#[cfg(test)]
mod test {
    use super::*;

    fn fixture_dir() -> PathBuf {
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data/pmsms_fixture")
    }

    #[test]
    fn mmappet_and_parquet_precursors_agree() {
        let dir = fixture_dir();
        let pmsms_dir = dir.join("pmsms.mmappet");

        let via_parquet = parse(&pmsms_dir, &dir.join("precursors.parquet"), 0)
            .expect("parse via parquet precursors");
        let via_mmappet = parse(&pmsms_dir, &dir.join("precursors.mmappet"), 0)
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
    }

    #[test]
    fn unrecognized_precursors_format_errors() {
        let dir = fixture_dir();
        let bogus = dir.join("precursors.parquet").with_extension("txt");
        let err = parse(&dir.join("pmsms.mmappet"), &bogus, 0).unwrap_err();
        assert!(matches!(err, PmsmsError::UnrecognizedPrecursorsFormat(_)));
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
            static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
            let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            let dir = std::env::temp_dir().join(format!(
                "sage_pmsms_test_{}_{name}_{n}.mmappet",
                std::process::id()
            ));
            let _ = std::fs::remove_dir_all(&dir);
            std::fs::create_dir_all(&dir).expect("create temp mmappet dir");
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

    fn mk_pmsms_fragments_dir(name: &str) -> PathBuf {
        // 2 precursors, 2 fragments each, values are arbitrary filler —
        // these tests only care about `isolation_window`.
        MmappetDirBuilder::new(name)
            .column("float32", "mz", &[100.0f32, 200.0, 300.0, 400.0])
            .column("uint32", "intensity", &[10u32, 20, 30, 40])
            .finish()
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
        let pmsms_dir = mk_pmsms_fragments_dir("pmsms_with_ppm");

        let spectra = parse(&pmsms_dir, &precursors_dir, 0).expect("parse with ppm columns");
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

        let _ = std::fs::remove_dir_all(&precursors_dir);
        let _ = std::fs::remove_dir_all(&pmsms_dir);
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
        let pmsms_dir = mk_pmsms_fragments_dir("pmsms_without_ppm");

        let spectra = parse(&pmsms_dir, &precursors_dir, 0).expect("parse without ppm columns");
        assert_eq!(spectra.len(), 2);
        for spectrum in &spectra {
            assert!(spectrum.precursors[0].isolation_window.is_none());
        }

        let _ = std::fs::remove_dir_all(&precursors_dir);
        let _ = std::fs::remove_dir_all(&pmsms_dir);
    }
}
