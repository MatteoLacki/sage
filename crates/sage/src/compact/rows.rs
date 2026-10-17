//! [`CompactPeptides`]' rows as columns: 12 B per row once a table is built, against the
//! 16 B [`CompactPeptide`] the builds work in (necromerge2 `plans/TODOS.md`, compact
//! table v2; `docs/ai/compact_peptides.md`).

use super::{CompactPeptide, VarSite, MAX_VARIABLE_MODS};
use rayon::prelude::*;

/// Rows per mass chunk.
const CHUNK: usize = 1024;

/// Rows sorted by precursor mass, one column per field. Row `i` decodes back to the
/// identical [`CompactPeptide`].
#[derive(Debug, Default)]
pub struct CompactRows {
    masses: ChunkedMasses,
    starts: Vec<u32>,
    lengths: Vec<u8>,
    flags: Vec<u8>,
    var_sites: Vec<[VarSite; MAX_VARIABLE_MODS]>,
}

impl CompactRows {
    /// Bytes a row takes in these columns (the chunk table adds ~8 B per 1,024 rows).
    pub const BYTES_PER_ROW: usize = 2 + 4 + 1 + 1 + 2 * MAX_VARIABLE_MODS;

    /// Columns of `rows`, which must be sorted by mass.
    pub fn from_rows(rows: Vec<CompactPeptide>) -> CompactRows {
        let masses = ChunkedMasses::new(&rows);
        let starts = rows.par_iter().map(|row| row.start).collect();
        let lengths = rows.par_iter().map(|row| row.length).collect();
        let flags = rows.par_iter().map(|row| row.flags).collect();
        let var_sites = rows.par_iter().map(|row| row.var_sites).collect();
        super::drop_in_parallel(rows);
        CompactRows {
            masses,
            starts,
            lengths,
            flags,
            var_sites,
        }
    }

    pub fn len(&self) -> usize {
        self.starts.len()
    }

    pub fn is_empty(&self) -> bool {
        self.starts.is_empty()
    }

    /// Row `row`, decoded.
    pub fn get(&self, row: usize) -> CompactPeptide {
        CompactPeptide {
            monoisotopic: self.masses.get(row),
            start: self.starts[row],
            length: self.lengths[row],
            flags: self.flags[row],
            var_sites: self.var_sites[row],
        }
    }

    /// Row `row`'s precursor mass.
    pub fn mass(&self, row: usize) -> f32 {
        self.masses.get(row)
    }

    /// Row `row`'s peptide length.
    pub fn length(&self, row: usize) -> usize {
        self.lengths[row] as usize
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
            + self.starts.len() * 4
            + self.lengths.len()
            + self.flags.len()
            + self.var_sites.len() * std::mem::size_of::<[VarSite; MAX_VARIABLE_MODS]>()
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
    residuals: Vec<u16>,
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
            residuals,
            wide,
        }
    }

    fn get(&self, row: usize) -> f32 {
        let (base, wide_start) = self.chunks[row / CHUNK];
        match wide_start {
            u32::MAX => f32::from_bits(base + self.residuals[row] as u32),
            _ => f32::from_bits(self.wide[wide_start as usize + row % CHUNK]),
        }
    }

    fn bytes(&self) -> usize {
        self.chunks.len() * 8 + self.residuals.len() * 2 + self.wide.len() * 4
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
}
