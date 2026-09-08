//! Cheap integrity check for the positional (`peptide_row`-indexed, not
//! `sequence`-keyed) `--predicted-rt`/`--predicted-iim`/
//! `--predicted-fragment-intensity-index` files `git/featureprediction`
//! exports from one `dumped_peptides` parquet -- see
//! `docs/ai/dumped_peptides_positional_predictions.md`.
//!
//! Row `i` of those files is only valid against a digest whose `PeptideIx`
//! assignment is byte-identical to the `dumped_peptides` node they were
//! built from (guaranteed when both share the same fasta + `database`
//! config, which necroflow's node-identity hashing already enforces --
//! see `git/sage`'s `docs/ai/dump_peptides.md`). This fingerprint is the
//! runtime check that guarantee actually held for *this* run, instead of
//! trusting it silently: a SHA-256 over each peptide's `monoisotopic` mass
//! (row order, little-endian `f32` bytes) is cheap to compute on both
//! sides -- no `Peptide::to_string()` (the expensive, unimod-lookup-per-
//! residue formatting this whole positional scheme exists to avoid) is
//! needed, `monoisotopic` is a plain field on both the Rust `Peptide` and
//! the `dumped_peptides` parquet's own `monoisotopic` column.

use sha2::{Digest, Sha256};

/// SHA-256 hex digest over `masses`' little-endian `f32` bytes, in order.
pub fn fingerprint_monoisotopic(masses: impl Iterator<Item = f32>) -> String {
    let mut hasher = Sha256::new();
    for m in masses {
        hasher.update(m.to_le_bytes());
    }
    format!("{:x}", hasher.finalize())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_and_order_sensitive() {
        let a = fingerprint_monoisotopic([1.0f32, 2.0, 3.0].into_iter());
        let b = fingerprint_monoisotopic([1.0f32, 2.0, 3.0].into_iter());
        let c = fingerprint_monoisotopic([3.0f32, 2.0, 1.0].into_iter());
        assert_eq!(a, b);
        assert_ne!(a, c);
    }
}
