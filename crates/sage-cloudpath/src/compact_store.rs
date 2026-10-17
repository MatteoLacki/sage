//! A [`CompactPeptides`] table stored once per FASTA and digest settings, and read
//! back in place: `dump_peptides --compact <dir>` writes it, `sage
//! --compact-peptides <dir>` maps it instead of building the table (see
//! `docs/ai/compact_peptides.md`).
//!
//! `<dir>` is an mmappet dataset over the rows (`schema.txt`, `0.bin`..`4.bin`:
//! the [`CompactRows`](sage_core::compact::CompactRows) columns, native byte order)
//! with three sidecars: `mass_chunks.bin` (per chunk of rows, the first mass's bit
//! pattern and where its full patterns start in `wide_masses.bin`, two `u32`),
//! `wide_masses.bin` (`u32` bit patterns) and `meta.json` ([`Meta`]). The proteome is
//! not stored: the reader rebuilds it from its own FASTA and checks it by hash.

use crate::dumped_peptides_fingerprint::fingerprint_monoisotopic;
use memmap2::Mmap;
use sage_core::compact::{CompactPeptides, Proteome, SharedBytes, StoredRows, ROWS_PER_MASS_CHUNK};
use sage_core::database::{EnzymeBuilder, Parameters};
use sage_core::fasta::Fasta;
use sage_core::modification::ModificationSpecificity;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;

#[cfg(not(target_endian = "little"))]
compile_error!("stored compact tables are written in little-endian byte order");

const FORMAT: &str = "sage_compact_peptides";
const VERSION: u32 = 1;
/// `schema.txt`: numpy dtype and name of each row column, in file order.
const ROW_COLUMNS: [(&str, &str); 5] = [
    ("uint16", "mass_residual"),
    ("uint32", "start"),
    ("uint8", "length"),
    ("uint8", "flags"),
    ("uint32", "var_sites"),
];
const MASS_CHUNKS: &str = "mass_chunks.bin";
const WIDE_MASSES: &str = "wide_masses.bin";
const META: &str = "meta.json";

#[derive(Debug, thiserror::Error)]
pub enum CompactStoreError {
    #[error("io error: {0}")]
    Io(#[from] std::io::Error),
    #[error("json error: {0}")]
    Json(#[from] serde_json::Error),
    #[error("{0}")]
    Table(String),
    #[error("stored compact peptides do not match this run: {0}")]
    Mismatch(String),
}

/// What a stored table records about itself, and what a reader checks.
#[derive(Serialize, Deserialize)]
struct Meta {
    format: String,
    version: u32,
    rows: usize,
    rows_per_mass_chunk: usize,
    /// `dumped_peptides_sha256` of the rows' masses, as positional predictions carry.
    dumped_peptides_sha256: String,
    proteome_sha256: String,
    /// The distinct variable-mod masses rows index, as `f32` bit patterns.
    variable_mod_mass_bits: Vec<u32>,
    /// The `database` settings the digest reads (see [`DigestSettings`]).
    digest_settings: serde_json::Value,
}

/// The `database` settings a digest depends on: the keys the pipeline's
/// `dump_peptides_config.json` holds, FASTA aside.
#[derive(Serialize)]
struct DigestSettings<'a> {
    decoy_tag: &'a str,
    enzyme: &'a EnzymeBuilder,
    generate_decoys: bool,
    max_variable_mods: usize,
    peptide_max_mass: f32,
    peptide_min_mass: f32,
    static_mods: &'a HashMap<ModificationSpecificity, f32>,
    variable_mods: &'a HashMap<ModificationSpecificity, Vec<f32>>,
}

/// `parameters`' digest settings as JSON, through text so both sides of a
/// comparison print `f32`s the same way.
fn digest_settings(parameters: &Parameters) -> Result<serde_json::Value, serde_json::Error> {
    let text = serde_json::to_string(&DigestSettings {
        decoy_tag: &parameters.decoy_tag,
        enzyme: &parameters.enzyme,
        generate_decoys: parameters.generate_decoys,
        max_variable_mods: parameters.max_variable_mods,
        peptide_max_mass: parameters.peptide_max_mass,
        peptide_min_mass: parameters.peptide_min_mass,
        static_mods: &parameters.static_mods,
        variable_mods: &parameters.variable_mods,
    })?;
    serde_json::from_str(&text)
}

/// SHA-256 over the proteome's protein boundaries, accessions and residues.
fn proteome_sha256(proteome: &Proteome) -> String {
    let mut hasher = Sha256::new();
    hasher.update((proteome.protein_count() as u64).to_le_bytes());
    for start in proteome.protein_starts() {
        hasher.update(start.to_le_bytes());
    }
    for accession in proteome.accessions() {
        hasher.update((accession.len() as u64).to_le_bytes());
        hasher.update(accession.as_bytes());
    }
    hasher.update(proteome.residues());
    format!("{:x}", hasher.finalize())
}

/// Store `table` (built with `parameters`) in `dir`. Returns the rows'
/// `dumped_peptides_sha256`.
pub fn write(
    table: &CompactPeptides,
    parameters: &Parameters,
    dir: &Path,
) -> Result<String, CompactStoreError> {
    std::fs::create_dir_all(dir)?;
    let columns = table.rows.columns();
    let schema: Vec<String> = ROW_COLUMNS
        .iter()
        .map(|(dtype, name)| format!("{dtype} {name}"))
        .collect();
    std::fs::write(dir.join("schema.txt"), schema.join("\n"))?;
    let row_bytes = [
        columns.mass_residuals,
        columns.starts,
        columns.lengths,
        columns.flags,
        columns.var_sites,
    ];
    for (index, bytes) in row_bytes.iter().enumerate() {
        std::fs::write(dir.join(format!("{index}.bin")), bytes)?;
    }
    let mass_chunks: Vec<u8> = columns
        .mass_chunks
        .iter()
        .flat_map(|&(base, wide_start)| [base.to_le_bytes(), wide_start.to_le_bytes()])
        .flatten()
        .collect();
    std::fs::write(dir.join(MASS_CHUNKS), mass_chunks)?;
    let wide_masses: Vec<u8> = columns
        .wide_masses
        .iter()
        .flat_map(|bits| bits.to_le_bytes())
        .collect();
    std::fs::write(dir.join(WIDE_MASSES), wide_masses)?;

    let fingerprint = fingerprint_monoisotopic(table.rows.masses());
    let meta = Meta {
        format: FORMAT.into(),
        version: VERSION,
        rows: columns.rows,
        rows_per_mass_chunk: ROWS_PER_MASS_CHUNK,
        dumped_peptides_sha256: fingerprint.clone(),
        proteome_sha256: proteome_sha256(&table.proteome),
        variable_mod_mass_bits: table
            .variable_mod_masses()
            .iter()
            .map(|mass| mass.to_bits())
            .collect(),
        digest_settings: digest_settings(parameters)?,
    };
    std::fs::write(dir.join(META), serde_json::to_vec_pretty(&meta)?)?;
    Ok(fingerprint)
}

/// The table stored in `dir`, its rows read in place, after checking that it was
/// built from `fasta` with `parameters`' digest settings.
pub fn open(
    dir: &Path,
    parameters: &Parameters,
    fasta: &Fasta,
) -> Result<CompactPeptides, CompactStoreError> {
    let meta: Meta = serde_json::from_slice(&std::fs::read(dir.join(META))?)?;
    if meta.format != FORMAT || meta.version != VERSION {
        return Err(CompactStoreError::Table(format!(
            "{} holds {} version {}, expected {FORMAT} version {VERSION}",
            dir.display(),
            meta.format,
            meta.version
        )));
    }
    if meta.rows_per_mass_chunk != ROWS_PER_MASS_CHUNK {
        return Err(CompactStoreError::Table(format!(
            "stored with {} rows per mass chunk, this build reads {ROWS_PER_MASS_CHUNK}",
            meta.rows_per_mass_chunk
        )));
    }
    let settings = digest_settings(parameters)?;
    if settings != meta.digest_settings {
        return Err(CompactStoreError::Mismatch(format!(
            "digest settings {} were stored, this run has {}",
            meta.digest_settings, settings
        )));
    }

    let mut table = CompactPeptides::new(parameters, fasta)
        .map_err(|e| CompactStoreError::Table(e.to_string()))?;
    let mod_bits: Vec<u32> = table
        .variable_mod_masses()
        .iter()
        .map(|mass| mass.to_bits())
        .collect();
    if mod_bits != meta.variable_mod_mass_bits {
        return Err(CompactStoreError::Mismatch(
            "the variable-mod masses differ".into(),
        ));
    }
    if proteome_sha256(&table.proteome) != meta.proteome_sha256 {
        return Err(CompactStoreError::Mismatch(
            "the FASTA's proteome differs from the one the table was built from".into(),
        ));
    }

    let schema = std::fs::read_to_string(dir.join("schema.txt"))?;
    let stored: Vec<(&str, &str)> = schema
        .lines()
        .filter(|line| !line.trim().is_empty())
        .filter_map(|line| line.split_once(' '))
        .collect();
    if stored != ROW_COLUMNS {
        return Err(CompactStoreError::Table(format!(
            "row columns {stored:?}, expected {ROW_COLUMNS:?}"
        )));
    }
    let mass_chunks = std::fs::read(dir.join(MASS_CHUNKS))?
        .chunks_exact(8)
        .map(|pair| {
            let base = u32::from_le_bytes(pair[..4].try_into().unwrap());
            let wide_start = u32::from_le_bytes(pair[4..].try_into().unwrap());
            (base, wide_start)
        })
        .collect();
    let wide_masses = std::fs::read(dir.join(WIDE_MASSES))?
        .chunks_exact(4)
        .map(|bits| u32::from_le_bytes(bits.try_into().unwrap()))
        .collect();
    table.rows = sage_core::compact::CompactRows::from_stored(StoredRows {
        rows: meta.rows,
        mass_chunks,
        wide_masses,
        mass_residuals: map_column(dir, 0)?,
        starts: map_column(dir, 1)?,
        lengths: map_column(dir, 2)?,
        flags: map_column(dir, 3)?,
        var_sites: map_column(dir, 4)?,
    })
    .map_err(CompactStoreError::Table)?;
    table.stored_fingerprint = Some(meta.dumped_peptides_sha256);
    Ok(table)
}

/// Row column `index`, mapped read-only (an empty file is an empty column).
fn map_column(dir: &Path, index: usize) -> Result<SharedBytes, CompactStoreError> {
    let file = std::fs::File::open(dir.join(format!("{index}.bin")))?;
    if file.metadata()?.len() == 0 {
        return Ok(Arc::new(Vec::<u8>::new()));
    }
    // SAFETY: the table's files are written once and never modified while mapped.
    let map = unsafe { Mmap::map(&file)? };
    Ok(Arc::new(map))
}

#[cfg(test)]
mod tests {
    use super::*;
    use sage_core::database::Builder;
    use std::path::PathBuf;

    const FASTA: &str = ">sp|P1|ONE\nMKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPF\n\
                         >sp|P2|TWO\nMDLSALRVEEVQNVINAMQKILECPICLELIKEPVSTKCDHIFCKFCMLKLLNQKKGPSQC\n\
                         >sp|P3|THREE\nMKWVTFISLLLLFSSAYSRGVFRR\n";

    fn temp_dir(name: &str) -> PathBuf {
        static COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
        let n = COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        let dir = std::env::temp_dir().join(format!(
            "sage_compact_store_test_{}_{name}_{n}",
            std::process::id()
        ));
        let _ = std::fs::remove_dir_all(&dir);
        dir
    }

    fn parameters(config: serde_json::Value) -> Parameters {
        let mut builder: Builder = serde_json::from_value(config).unwrap();
        builder.fasta = Some("unused".into());
        builder.make_parameters()
    }

    fn nonspecific() -> serde_json::Value {
        serde_json::json!({
            "enzyme": {"cleave_at": "", "min_len": 7, "max_len": 11},
            "peptide_min_mass": 500.0,
            "peptide_max_mass": 2000.0,
            "static_mods": {},
            "variable_mods": {"M": [15.9949], "[": [42.010567]},
            "max_variable_mods": 2,
            "decoy_tag": "rev_",
            "generate_decoys": true
        })
    }

    fn tryptic() -> serde_json::Value {
        serde_json::json!({
            "enzyme": {"missed_cleavages": 2, "min_len": 5, "max_len": 30},
            "static_mods": {"C": 57.021464},
            "variable_mods": {"M": [15.9949]},
            "decoy_tag": "rev_",
            "generate_decoys": true
        })
    }

    fn assert_round_trip(config: serde_json::Value) {
        let parameters = parameters(config);
        let fasta = Fasta::parse(FASTA.into(), "rev_", true);
        let built = CompactPeptides::build(&parameters, &fasta).unwrap();
        assert!(!built.rows.is_empty());
        let dir = temp_dir("round_trip");
        let fingerprint = write(&built, &parameters, &dir).unwrap();

        let opened = open(&dir, &parameters, &fasta).unwrap();
        assert!(opened.rows.is_stored());
        assert_eq!(
            opened.stored_fingerprint.as_deref(),
            Some(fingerprint.as_str())
        );
        assert_eq!(
            fingerprint,
            fingerprint_monoisotopic(opened.rows.masses()),
            "the stored fingerprint is the stored rows'"
        );
        assert_eq!(opened.rows.len(), built.rows.len());
        for row in 0..built.rows.len() {
            assert_eq!(opened.rows.get(row), built.rows.get(row), "row {row}");
            assert_eq!(
                opened.materialize(&opened.rows.get(row)).to_string(),
                built.materialize(&built.rows.get(row)).to_string()
            );
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn nonspecific_table_round_trips() {
        assert_round_trip(nonspecific());
    }

    #[test]
    fn tryptic_table_round_trips() {
        assert_round_trip(tryptic());
    }

    fn open_error(dir: &Path, parameters: &Parameters, fasta: &Fasta) -> CompactStoreError {
        match open(dir, parameters, fasta) {
            Ok(_) => panic!("opened a table stored for other settings"),
            Err(err) => err,
        }
    }

    #[test]
    fn rejects_other_settings_and_proteomes() {
        let stored_with = parameters(nonspecific());
        let fasta = Fasta::parse(FASTA.into(), "rev_", true);
        let built = CompactPeptides::build(&stored_with, &fasta).unwrap();
        let dir = temp_dir("rejects");
        write(&built, &stored_with, &dir).unwrap();

        let mut other_mods = nonspecific();
        other_mods["variable_mods"] = serde_json::json!({"M": [15.9949]});
        let err = open_error(&dir, &parameters(other_mods), &fasta);
        assert!(matches!(err, CompactStoreError::Mismatch(_)), "{err}");

        let mut other_lengths = nonspecific();
        other_lengths["enzyme"]["max_len"] = serde_json::json!(12);
        let err = open_error(&dir, &parameters(other_lengths), &fasta);
        assert!(matches!(err, CompactStoreError::Mismatch(_)), "{err}");

        let other_fasta = Fasta::parse(FASTA.replace("QYLQQ", "QYLQE"), "rev_", true);
        let err = open_error(&dir, &stored_with, &other_fasta);
        assert!(matches!(err, CompactStoreError::Mismatch(_)), "{err}");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
