//! [`CompactPeptides`]' rows as columns: 12 B per row once a table is built, against the
//! 16 B [`CompactPeptide`] the builds work in (necromerge2 `plans/TODOS.md`, compact
//! table v2; `docs/ai/compact_peptides.md`).

use super::{CompactPeptide, VarSite, MAX_VARIABLE_MODS};
use rayon::prelude::*;
use std::sync::Arc;

/// Rows per mass chunk.
pub const CHUNK: usize = 1024;

type VarSites = [VarSite; MAX_VARIABLE_MODS];

/// Bytes kept alive for as long as a column is read from them, e.g. a memory-mapped
/// file.
pub type SharedBytes = Arc<dyn AsRef<[u8]> + Send + Sync>;

/// Element types a column can be read from stored bytes as: every bit pattern of
/// their size is a valid value.
trait StoredElement: Copy + Send + Sync + 'static {}
impl StoredElement for u8 {}
impl StoredElement for u16 {}
impl StoredElement for u32 {}
impl StoredElement for VarSites {}

/// One column: built in memory, or read in place from stored bytes. Either way it is
/// read through the same pointer and length, so a stored column costs no more per
/// access than an owned one.
struct Column<T: StoredElement> {
    values: *const T,
    len: usize,
    storage: Storage<T>,
}

/// What keeps a [`Column`]'s values alive (never read: the column reads through its
/// pointer).
enum Storage<T> {
    Owned { _values: Vec<T> },
    Stored { _bytes: SharedBytes },
}

// SAFETY: a column's values are never written after construction, and `storage`
// (itself `Send + Sync`) keeps them alive and in place.
unsafe impl<T: StoredElement> Send for Column<T> {}
unsafe impl<T: StoredElement> Sync for Column<T> {}

impl<T: StoredElement> Column<T> {
    fn owned(values: Vec<T>) -> Self {
        Self {
            values: values.as_ptr(),
            len: values.len(),
            storage: Storage::Owned { _values: values },
        }
    }

    /// `bytes` read as `rows` values of `T`, native byte order.
    fn stored(bytes: SharedBytes, rows: usize, name: &str) -> Result<Self, String> {
        let raw: &[u8] = (*bytes).as_ref();
        let expected = rows * std::mem::size_of::<T>();
        if raw.len() != expected {
            return Err(format!(
                "column {name} holds {} bytes, expected {expected} for {rows} rows",
                raw.len()
            ));
        }
        let values = match rows {
            0 => std::ptr::NonNull::<T>::dangling().as_ptr() as *const T,
            _ => raw.as_ptr() as *const T,
        };
        if values as usize % std::mem::align_of::<T>() != 0 {
            return Err(format!("column {name} is not aligned for its element type"));
        }
        Ok(Self {
            values,
            len: rows,
            storage: Storage::Stored { _bytes: bytes },
        })
    }

    fn as_slice(&self) -> &[T] {
        // SAFETY: `values` points at `len` initialised values of `T` (checked for
        // length and alignment when stored), kept alive by `storage` for as long as
        // `self`.
        unsafe { std::slice::from_raw_parts(self.values, self.len) }
    }

    fn is_stored(&self) -> bool {
        matches!(self.storage, Storage::Stored { .. })
    }

    fn bytes(&self) -> &[u8] {
        let values = self.as_slice();
        // SAFETY: plain values without padding viewed as their bytes.
        unsafe {
            std::slice::from_raw_parts(values.as_ptr() as *const u8, std::mem::size_of_val(values))
        }
    }
}

impl<T: StoredElement> Default for Column<T> {
    fn default() -> Self {
        Column::owned(Vec::new())
    }
}

impl<T: StoredElement> std::fmt::Debug for Column<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let kind = match self.storage {
            Storage::Owned { .. } => "Owned",
            Storage::Stored { .. } => "Stored",
        };
        write!(f, "{kind}({} values)", self.len)
    }
}

/// [`CompactRows`]' columns as bytes, to store; [`StoredRows`] reads them back.
pub struct RowColumns<'a> {
    pub rows: usize,
    /// Per [`CHUNK`] rows: the first row's mass bit pattern, and where the chunk's
    /// full bit patterns start in `wide_masses` (`u32::MAX` when it has none).
    pub mass_chunks: &'a [(u32, u32)],
    pub wide_masses: &'a [u32],
    pub mass_residuals: &'a [u8],
    pub starts: &'a [u8],
    pub lengths: &'a [u8],
    pub flags: &'a [u8],
    pub var_sites: &'a [u8],
}

/// Stored columns, as [`RowColumns`] wrote them, for [`CompactRows::from_stored`].
pub struct StoredRows {
    pub rows: usize,
    pub mass_chunks: Vec<(u32, u32)>,
    pub wide_masses: Vec<u32>,
    pub mass_residuals: SharedBytes,
    pub starts: SharedBytes,
    pub lengths: SharedBytes,
    pub flags: SharedBytes,
    pub var_sites: SharedBytes,
}

/// Rows sorted by precursor mass, one column per field. Row `i` decodes back to the
/// identical [`CompactPeptide`].
#[derive(Debug, Default)]
pub struct CompactRows {
    masses: ChunkedMasses,
    starts: Column<u32>,
    lengths: Column<u8>,
    flags: Column<u8>,
    var_sites: Column<VarSites>,
}

impl CompactRows {
    /// Bytes a row takes in these columns (the chunk table adds ~8 B per 1,024 rows).
    pub const BYTES_PER_ROW: usize = 2 + 4 + 1 + 1 + 2 * MAX_VARIABLE_MODS;

    /// Columns of `rows`, which must be sorted by mass.
    pub fn from_rows(rows: Vec<CompactPeptide>) -> CompactRows {
        let masses = ChunkedMasses::new(&rows);
        let starts = Column::owned(rows.par_iter().map(|row| row.start).collect());
        let lengths = Column::owned(rows.par_iter().map(|row| row.length).collect());
        let flags = Column::owned(rows.par_iter().map(|row| row.flags).collect());
        let var_sites = Column::owned(rows.par_iter().map(|row| row.var_sites).collect());
        super::drop_in_parallel(rows);
        CompactRows {
            masses,
            starts,
            lengths,
            flags,
            var_sites,
        }
    }

    /// Rows read in place from stored columns, checked for size and consistency.
    pub fn from_stored(stored: StoredRows) -> Result<CompactRows, String> {
        let rows = stored.rows;
        let masses = ChunkedMasses::stored(
            rows,
            stored.mass_chunks,
            stored.wide_masses,
            Column::stored(stored.mass_residuals, rows, "mass_residual")?,
        )?;
        Ok(CompactRows {
            masses,
            starts: Column::stored(stored.starts, rows, "start")?,
            lengths: Column::stored(stored.lengths, rows, "length")?,
            flags: Column::stored(stored.flags, rows, "flags")?,
            var_sites: Column::stored(stored.var_sites, rows, "var_sites")?,
        })
    }

    /// The columns as bytes, in native byte order.
    pub fn columns(&self) -> RowColumns<'_> {
        RowColumns {
            rows: self.len(),
            mass_chunks: &self.masses.chunks,
            wide_masses: &self.masses.wide,
            mass_residuals: self.masses.residuals.bytes(),
            starts: self.starts.bytes(),
            lengths: self.lengths.bytes(),
            flags: self.flags.bytes(),
            var_sites: self.var_sites.bytes(),
        }
    }

    pub fn len(&self) -> usize {
        self.starts.as_slice().len()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Row `row`, decoded.
    pub fn get(&self, row: usize) -> CompactPeptide {
        CompactPeptide {
            monoisotopic: self.masses.get(row),
            start: self.starts.as_slice()[row],
            length: self.lengths.as_slice()[row],
            flags: self.flags.as_slice()[row],
            var_sites: self.var_sites.as_slice()[row],
        }
    }

    /// Row `row`'s precursor mass.
    pub fn mass(&self, row: usize) -> f32 {
        self.masses.get(row)
    }

    /// Row `row`'s peptide length.
    pub fn length(&self, row: usize) -> usize {
        self.lengths.as_slice()[row] as usize
    }

    pub fn iter(&self) -> impl Iterator<Item = CompactPeptide> + '_ {
        (0..self.len()).map(|row| self.get(row))
    }

    pub fn masses(&self) -> impl Iterator<Item = f32> + '_ {
        (0..self.len()).map(|row| self.masses.get(row))
    }

    /// Bytes the columns hold.
    pub fn bytes(&self) -> usize {
        self.masses.bytes()
            + self.starts.bytes().len()
            + self.lengths.bytes().len()
            + self.flags.bytes().len()
            + self.var_sites.bytes().len()
    }

    /// Whether the columns are read in place from stored bytes.
    pub fn is_stored(&self) -> bool {
        self.starts.is_stored()
    }
}

/// Masses sorted ascending, as 16-bit offsets from their chunk's first f32 bit pattern:
/// positive floats order as their bits, so `from_bits(base + residual)` is the exact
/// mass. A chunk spanning more than 2^16 bit patterns (only in the sparse mass tails)
/// keeps its rows' full bit patterns instead.
#[derive(Debug, Default)]
struct ChunkedMasses {
    /// Per chunk: its first row's bit pattern, and where its full patterns start in
    /// `wide` (`u32::MAX` for a chunk that fits in residuals).
    chunks: Vec<(u32, u32)>,
    residuals: Column<u16>,
    wide: Vec<u32>,
}

impl ChunkedMasses {
    fn new(rows: &[CompactPeptide]) -> ChunkedMasses {
        let bits_of = |row: &CompactPeptide| {
            assert!(
                row.monoisotopic >= 0.0 && row.monoisotopic.is_finite(),
                "row mass {} cannot be stored as a residual",
                row.monoisotopic
            );
            row.monoisotopic.to_bits()
        };
        let mut chunks = Vec::with_capacity((rows.len() + CHUNK - 1) / CHUNK);
        let mut wide = Vec::new();
        for chunk in rows.chunks(CHUNK) {
            let base = bits_of(&chunk[0]);
            let last = bits_of(&chunk[chunk.len() - 1]);
            match last - base <= u16::MAX as u32 {
                true => chunks.push((base, u32::MAX)),
                false => {
                    chunks.push((base, wide.len() as u32));
                    wide.extend(chunk.iter().map(bits_of));
                }
            }
        }
        let residuals = rows
            .par_iter()
            .enumerate()
            .map(|(row, peptide)| {
                let (base, wide_start) = chunks[row / CHUNK];
                match wide_start {
                    u32::MAX => (peptide.monoisotopic.to_bits() - base) as u16,
                    _ => 0,
                }
            })
            .collect();
        ChunkedMasses {
            chunks,
            residuals: Column::owned(residuals),
            wide,
        }
    }

    /// Stored chunks over `rows` rows, checked so every row decodes inside them.
    fn stored(
        rows: usize,
        chunks: Vec<(u32, u32)>,
        wide: Vec<u32>,
        residuals: Column<u16>,
    ) -> Result<ChunkedMasses, String> {
        if chunks.len() != (rows + CHUNK - 1) / CHUNK {
            return Err(format!(
                "{} mass chunks stored for {rows} rows, expected {}",
                chunks.len(),
                (rows + CHUNK - 1) / CHUNK
            ));
        }
        for (chunk, &(_, wide_start)) in chunks.iter().enumerate() {
            let chunk_rows = CHUNK.min(rows - chunk * CHUNK);
            if wide_start != u32::MAX && wide_start as usize + chunk_rows > wide.len() {
                return Err(format!("mass chunk {chunk} points past the wide masses"));
            }
        }
        Ok(ChunkedMasses {
            chunks,
            residuals,
            wide,
        })
    }

    fn get(&self, row: usize) -> f32 {
        let (base, wide_start) = self.chunks[row / CHUNK];
        match wide_start {
            u32::MAX => f32::from_bits(base + self.residuals.as_slice()[row] as u32),
            _ => f32::from_bits(self.wide[wide_start as usize + row % CHUNK]),
        }
    }

    fn bytes(&self) -> usize {
        self.chunks.len() * 8 + self.residuals.bytes().len() + self.wide.len() * 4
    }
}

#[cfg(test)]
mod test {
    use super::*;

    fn row(monoisotopic: f32, start: u32) -> CompactPeptide {
        CompactPeptide {
            monoisotopic,
            start,
            length: (start % 20) as u8 + 5,
            flags: (start % 7) as u8,
            var_sites: [
                VarSite::new(crate::peptide::Site::Sequence(start % 5), (start % 3) as u8),
                VarSite::UNUSED,
            ],
        }
    }

    #[test]
    fn columns_decode_every_row_exactly_including_wide_chunks() {
        // Dense masses, then a sparse tail whose chunks span far more than 2^16 patterns.
        let mut masses: Vec<f32> = (0..5_000).map(|i| 800.0 + i as f32 * 1e-4).collect();
        masses.extend((0..3_000).map(|i| 2_000.0 + i as f32 * 7.5));
        masses.push(masses[masses.len() - 1]);
        let rows: Vec<CompactPeptide> = masses
            .iter()
            .enumerate()
            .map(|(i, &m)| row(m, i as u32))
            .collect();
        let columns = CompactRows::from_rows(rows.clone());
        assert_eq!(columns.len(), rows.len());
        assert!(
            !columns.masses.wide.is_empty(),
            "the tail needs wide chunks"
        );
        assert!(columns.masses.chunks.iter().any(|&(_, w)| w == u32::MAX));
        for (i, expected) in rows.iter().enumerate() {
            let decoded = columns.get(i);
            assert_eq!(decoded, *expected, "row {i}");
            assert_eq!(
                decoded.monoisotopic.to_bits(),
                expected.monoisotopic.to_bits()
            );
        }
        assert!(columns
            .masses()
            .zip(&masses)
            .all(|(a, b)| a.to_bits() == b.to_bits()));
        assert!(CompactRows::from_rows(Vec::new()).is_empty());
    }

    /// Bytes in a `u64` buffer, so every column type is aligned in them.
    struct AlignedBytes(Vec<u64>, usize);

    impl AsRef<[u8]> for AlignedBytes {
        fn as_ref(&self) -> &[u8] {
            // SAFETY: the buffer holds at least `self.1` initialised bytes.
            unsafe { std::slice::from_raw_parts(self.0.as_ptr() as *const u8, self.1) }
        }
    }

    fn shared(bytes: &[u8]) -> SharedBytes {
        let mut words = vec![0u64; (bytes.len() + 7) / 8];
        // SAFETY: `words` has room for `bytes.len()` bytes.
        unsafe {
            std::ptr::copy_nonoverlapping(
                bytes.as_ptr(),
                words.as_mut_ptr() as *mut u8,
                bytes.len(),
            )
        };
        Arc::new(AlignedBytes(words, bytes.len()))
    }

    fn stored_copy(columns: &CompactRows) -> Result<CompactRows, String> {
        let c = columns.columns();
        CompactRows::from_stored(StoredRows {
            rows: c.rows,
            mass_chunks: c.mass_chunks.to_vec(),
            wide_masses: c.wide_masses.to_vec(),
            mass_residuals: shared(c.mass_residuals),
            starts: shared(c.starts),
            lengths: shared(c.lengths),
            flags: shared(c.flags),
            var_sites: shared(c.var_sites),
        })
    }

    #[test]
    fn stored_columns_decode_every_row_exactly() {
        let mut masses: Vec<f32> = (0..5_000).map(|i| 800.0 + i as f32 * 1e-4).collect();
        masses.extend((0..3_000).map(|i| 2_000.0 + i as f32 * 7.5));
        let rows: Vec<CompactPeptide> = masses
            .iter()
            .enumerate()
            .map(|(i, &m)| row(m, i as u32))
            .collect();
        let columns = CompactRows::from_rows(rows.clone());
        let stored = stored_copy(&columns).unwrap();
        assert!(stored.is_stored() && !columns.is_stored());
        assert_eq!(stored.len(), rows.len());
        assert_eq!(stored.bytes(), columns.bytes());
        for (i, expected) in rows.iter().enumerate() {
            assert_eq!(stored.get(i), *expected, "row {i}");
            assert_eq!(stored.mass(i).to_bits(), expected.monoisotopic.to_bits());
        }
        assert!(stored_copy(&CompactRows::from_rows(Vec::new()))
            .unwrap()
            .is_empty());
        let empty: SharedBytes = Arc::new(Vec::<u8>::new());
        let empty_rows = CompactRows::from_stored(StoredRows {
            rows: 0,
            mass_chunks: Vec::new(),
            wide_masses: Vec::new(),
            mass_residuals: empty.clone(),
            starts: empty.clone(),
            lengths: empty.clone(),
            flags: empty.clone(),
            var_sites: empty,
        })
        .unwrap();
        assert!(empty_rows.is_empty() && empty_rows.columns().starts.is_empty());
    }

    #[test]
    fn stored_columns_are_checked() {
        let rows: Vec<CompactPeptide> = (0..3_000).map(|i| row(900.0 + i as f32, i)).collect();
        let columns = CompactRows::from_rows(rows);
        let c = columns.columns();
        let mut short_starts = c.starts.to_vec();
        short_starts.truncate(short_starts.len() - 4);
        let wrong_length = CompactRows::from_stored(StoredRows {
            rows: c.rows,
            mass_chunks: c.mass_chunks.to_vec(),
            wide_masses: c.wide_masses.to_vec(),
            mass_residuals: shared(c.mass_residuals),
            starts: shared(&short_starts),
            lengths: shared(c.lengths),
            flags: shared(c.flags),
            var_sites: shared(c.var_sites),
        });
        assert!(wrong_length.unwrap_err().contains("start"));
        let missing_chunk = CompactRows::from_stored(StoredRows {
            rows: c.rows,
            mass_chunks: c.mass_chunks[1..].to_vec(),
            wide_masses: c.wide_masses.to_vec(),
            mass_residuals: shared(c.mass_residuals),
            starts: shared(c.starts),
            lengths: shared(c.lengths),
            flags: shared(c.flags),
            var_sites: shared(c.var_sites),
        });
        assert!(missing_chunk.unwrap_err().contains("mass chunks"));
    }
}
