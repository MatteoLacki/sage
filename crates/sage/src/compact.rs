//! A peptide database stored as 16 B rows over one shared copy of the proteome,
//! in exactly the row order [`Parameters::digest`] produces, so row `i` is the
//! `PeptideIx(i)` a monolithic search would assign. A full [`Peptide`] is only
//! created on demand, by [`CompactPeptides::materialize`].

use crate::database::{IndexedDatabase, Parameters, ReportedPeptides};
use crate::enzyme::{Digest, EnzymeParameters, Position};
use crate::fasta::Fasta;
use crate::peptide::{Peptide, Site};
use crate::peptide_view::PeptideView;
use dashmap::DashMap;
use fnv::{FnvBuildHasher, FnvHashMap};
use rayon::prelude::*;
use std::borrow::Cow;
use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::sync::atomic::{AtomicU32, Ordering as AtomicOrdering};
use std::sync::Arc;

mod nonspecific;
mod rows;
mod view;

pub use rows::{CompactRows, RowColumns, SharedBytes, StoredRows, CHUNK as ROWS_PER_MASS_CHUNK};
pub use view::CompactView;

/// Residue indices share a `u8` with the two peptide-terminus codes.
const MAX_PEPTIDE_LENGTH: usize = VarSite::NTERM as usize;
const MAX_VARIABLE_MODS: usize = 2;
const MAX_MISSED_CLEAVAGES: u8 = 15;

#[derive(Debug, Clone, PartialEq)]
pub enum CompactError {
    DecoysNotGenerated,
    TooManyVariableMods(usize),
    TooManyVariableModMasses(usize),
    PeptideTooLong(usize),
    TooManyResidues,
    TooManyPeptides(usize),
    TooManyMissedCleavages(u8),
    UnrepresentableModification(String),
}

impl std::fmt::Display for CompactError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CompactError::DecoysNotGenerated => {
                write!(f, "compact peptides require `generate_decoys: true`")
            }
            CompactError::TooManyVariableMods(n) => write!(
                f,
                "compact peptides hold at most {MAX_VARIABLE_MODS} variable mods per peptide, got {n}"
            ),
            CompactError::TooManyVariableModMasses(n) => {
                write!(f, "compact peptides hold at most 254 distinct variable-mod masses, got {n}")
            }
            CompactError::PeptideTooLong(n) => write!(
                f,
                "compact peptides hold at most {MAX_PEPTIDE_LENGTH} residues, got max_len {n}"
            ),
            CompactError::TooManyResidues => {
                write!(f, "the FASTA's target residues do not fit a u32 offset")
            }
            CompactError::TooManyPeptides(n) => {
                write!(f, "{n} peptides do not fit a u32 peptide index")
            }
            CompactError::TooManyMissedCleavages(n) => write!(
                f,
                "compact peptides hold at most {MAX_MISSED_CLEAVAGES} missed cleavages, got {n}"
            ),
            CompactError::UnrepresentableModification(peptide) => {
                write!(f, "cannot express the modifications of {peptide} as variable-mod sites")
            }
        }
    }
}

impl std::error::Error for CompactError {}

/// Every target protein's residues, concatenated, so a peptide is an offset plus a length.
pub struct Proteome {
    residues: Vec<u8>,
    /// Protein `i` occupies `residues[protein_starts[i]..protein_starts[i + 1]]`.
    protein_starts: Vec<u32>,
    accessions: Vec<Arc<str>>,
}

impl Proteome {
    fn new(fasta: &Fasta) -> Result<Self, CompactError> {
        let mut residues = Vec::new();
        let mut protein_starts = vec![0];
        let mut accessions = Vec::with_capacity(fasta.targets.len());
        for (accession, sequence) in &fasta.targets {
            residues.extend_from_slice(sequence.as_bytes());
            let end = u32::try_from(residues.len()).map_err(|_| CompactError::TooManyResidues)?;
            protein_starts.push(end);
            accessions.push(accession.clone());
        }
        Ok(Self {
            residues,
            protein_starts,
            accessions,
        })
    }

    pub fn protein_count(&self) -> usize {
        self.accessions.len()
    }

    pub fn residue_count(&self) -> usize {
        self.residues.len()
    }

    /// Every target protein's residues, concatenated.
    pub fn residues(&self) -> &[u8] {
        &self.residues
    }

    /// Where each protein starts in [`Proteome::residues`], then where the last ends.
    pub fn protein_starts(&self) -> &[u32] {
        &self.protein_starts
    }

    pub fn accessions(&self) -> &[Arc<str>] {
        &self.accessions
    }

    fn protein_range(&self, protein: usize) -> std::ops::Range<usize> {
        self.protein_starts[protein] as usize..self.protein_starts[protein + 1] as usize
    }

    fn protein_sequence(&self, protein: usize) -> &str {
        std::str::from_utf8(&self.residues[self.protein_range(protein)])
            .expect("each protein's residues are a whole FASTA string")
    }

    fn window(&self, start: u32, length: u8) -> &[u8] {
        &self.residues[start as usize..start as usize + length as usize]
    }

    /// For each sequence, the proteome offset of its first occurrence inside a protein.
    fn first_window_starts(&self, sequences: &[Arc<[u8]>]) -> Vec<u32> {
        let distinct_sequences = AtomicU32::new(0);
        let slot_of_sequence: DashMap<&[u8], u32, FnvBuildHasher> =
            DashMap::with_capacity_and_hasher(sequences.len() / 2, FnvBuildHasher::default());
        let slots: Vec<u32> = sequences
            .par_iter()
            .map(|sequence| {
                *slot_of_sequence
                    .entry(sequence.as_ref())
                    .or_insert_with(|| distinct_sequences.fetch_add(1, AtomicOrdering::Relaxed))
            })
            .collect();
        let slot_of_sequence = slot_of_sequence.into_read_only();
        let lengths: BTreeSet<usize> = sequences.par_iter().map(|s| s.len()).collect();

        let first_starts: Vec<AtomicU32> = (0..distinct_sequences.into_inner())
            .map(|_| AtomicU32::new(u32::MAX))
            .collect();
        (0..self.protein_count()).into_par_iter().for_each(|protein| {
            let range = self.protein_range(protein);
            let residues = &self.residues[range.clone()];
            for &length in &lengths {
                for (offset, window) in residues.windows(length).enumerate() {
                    if let Some(&slot) = slot_of_sequence.get(window) {
                        first_starts[slot as usize]
                            .fetch_min((range.start + offset) as u32, AtomicOrdering::Relaxed);
                    }
                }
            }
        });

        slots
            .iter()
            .map(|&slot| {
                let start = first_starts[slot as usize].load(AtomicOrdering::Relaxed);
                assert_ne!(
                    start,
                    u32::MAX,
                    "every digested sequence is a window of some target protein"
                );
                start
            })
            .collect()
    }
}

/// One applied variable modification: a residue index or peptide terminus, and an
/// index into [`CompactPeptides`]' distinct variable-mod masses. Two bytes, any
/// value of either a valid one, so a stored column of them can be mapped as is.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
#[repr(C)]
struct VarSite {
    position: u8,
    mass_id: u8,
}

impl VarSite {
    const NTERM: u8 = u8::MAX - 1;
    const CTERM: u8 = u8::MAX;
    const UNUSED_MASS_ID: u8 = u8::MAX;
    const UNUSED: VarSite = VarSite {
        position: 0,
        mass_id: Self::UNUSED_MASS_ID,
    };

    fn new(site: Site, mass_id: u8) -> Self {
        let position = match site {
            Site::Nterm => Self::NTERM,
            Site::Cterm => Self::CTERM,
            Site::Sequence(index) => index as u8,
        };
        Self { position, mass_id }
    }

    fn site(self) -> Option<Site> {
        if self.mass_id == Self::UNUSED_MASS_ID {
            return None;
        }
        Some(match self.position {
            Self::NTERM => Site::Nterm,
            Self::CTERM => Site::Cterm,
            index => Site::Sequence(index as u32),
        })
    }
}

impl Default for VarSite {
    fn default() -> Self {
        Self::UNUSED
    }
}

/// One peptide in 16 B. Static mods, monoisotopic mass recomputation and the decoy
/// sequence are derived on demand; proteins are looked up only when asked for.
#[derive(Copy, Clone, Debug, PartialEq, Default)]
pub struct CompactPeptide {
    pub monoisotopic: f32,
    /// Proteome offset of the target sequence; a decoy reads it through [`Peptide::reverse`].
    start: u32,
    length: u8,
    /// [`CompactPeptide::DECOY`], [`CompactPeptide::SEMI_ENZYMATIC`], position, missed cleavages.
    flags: u8,
    var_sites: [VarSite; MAX_VARIABLE_MODS],
}

impl CompactPeptide {
    const DECOY: u8 = 1;
    const SEMI_ENZYMATIC: u8 = 1 << 1;
    const POSITION_SHIFT: u8 = 2;
    const POSITION_MASK: u8 = 0b11;
    const MISSED_CLEAVAGES_SHIFT: u8 = 4;

    fn new(
        monoisotopic: f32,
        start: u32,
        target: &Peptide,
    ) -> Result<CompactPeptide, CompactError> {
        if target.missed_cleavages > MAX_MISSED_CLEAVAGES {
            return Err(CompactError::TooManyMissedCleavages(target.missed_cleavages));
        }
        let mut flags = (Self::position_code(target.position) << Self::POSITION_SHIFT)
            | (target.missed_cleavages << Self::MISSED_CLEAVAGES_SHIFT);
        if target.semi_enzymatic {
            flags |= Self::SEMI_ENZYMATIC;
        }
        Ok(CompactPeptide {
            monoisotopic,
            start,
            length: target.sequence.len() as u8,
            flags,
            var_sites: [VarSite::UNUSED; MAX_VARIABLE_MODS],
        })
    }

    /// A non-specific digest's target row: no missed cleavages, not semi-enzymatic.
    fn nonspecific_target(
        monoisotopic: f32,
        start: u32,
        length: usize,
        position: Position,
        var_sites: [VarSite; MAX_VARIABLE_MODS],
    ) -> CompactPeptide {
        CompactPeptide {
            monoisotopic,
            start,
            length: length as u8,
            flags: Self::position_code(position) << Self::POSITION_SHIFT,
            var_sites,
        }
    }

    /// The decoy `Peptide::reverse` makes of this target row: same start, mass and sites.
    fn decoy_of(self) -> CompactPeptide {
        CompactPeptide {
            flags: self.flags | Self::DECOY,
            ..self
        }
    }

    /// The target row this decoy row reverses; a target row itself.
    fn target_of(self) -> CompactPeptide {
        CompactPeptide {
            flags: self.flags & !Self::DECOY,
            ..self
        }
    }

    fn position_code(position: Position) -> u8 {
        match position {
            Position::Nterm => 0,
            Position::Cterm => 1,
            Position::Full => 2,
            Position::Internal => 3,
        }
    }

    pub fn len(&self) -> usize {
        self.length as usize
    }

    pub fn is_decoy(&self) -> bool {
        self.flags & Self::DECOY != 0
    }

    fn is_semi_enzymatic(&self) -> bool {
        self.flags & Self::SEMI_ENZYMATIC != 0
    }

    /// Proteome offset of one occurrence of the target sequence.
    pub fn start(&self) -> u32 {
        self.start
    }

    pub fn position(&self) -> Position {
        match (self.flags >> Self::POSITION_SHIFT) & Self::POSITION_MASK {
            0 => Position::Nterm,
            1 => Position::Cterm,
            2 => Position::Full,
            _ => Position::Internal,
        }
    }

    fn missed_cleavages(&self) -> u8 {
        self.flags >> Self::MISSED_CLEAVAGES_SHIFT
    }
}

pub struct CompactPeptides {
    pub proteome: Proteome,
    /// `Parameters::digest`'s own order: row `i` is `PeptideIx(i)` of a monolithic search.
    pub rows: CompactRows,
    /// Distinct variable-mod masses; a [`VarSite`] stores an index into this.
    variable_mod_masses: Vec<f32>,
    parameters: Parameters,
    /// `parameters.static_mods` in the order `materialize` iterates them.
    static_mods: view::StaticMods,
    no_protein: Arc<str>,
    /// `dumped_peptides_sha256` of the rows' masses, recorded when the table was
    /// stored, so a run reading it back need not hash every mass again.
    pub stored_fingerprint: Option<String>,
}

impl CompactPeptides {
    /// An empty table over `fasta`'s proteome; [`CompactPeptides::build`] fills it.
    pub fn new(parameters: &Parameters, fasta: &Fasta) -> Result<Self, CompactError> {
        if !parameters.generate_decoys {
            return Err(CompactError::DecoysNotGenerated);
        }
        let variable_mods = parameters.variable_mod_list();
        if !variable_mods.is_empty() && parameters.max_variable_mods > MAX_VARIABLE_MODS {
            return Err(CompactError::TooManyVariableMods(parameters.max_variable_mods));
        }
        let enzyme: EnzymeParameters = parameters.enzyme.clone().into();
        if enzyme.max_len > MAX_PEPTIDE_LENGTH {
            return Err(CompactError::PeptideTooLong(enzyme.max_len));
        }

        let mut variable_mod_masses: Vec<f32> =
            variable_mods.into_iter().map(|(_, mass)| mass).collect();
        variable_mod_masses.sort_by(f32::total_cmp);
        variable_mod_masses.dedup_by(|a, b| a.to_bits() == b.to_bits());
        if variable_mod_masses.len() >= VarSite::UNUSED_MASS_ID as usize {
            return Err(CompactError::TooManyVariableModMasses(
                variable_mod_masses.len(),
            ));
        }

        let parameters = parameters.clone();
        Ok(Self {
            proteome: Proteome::new(fasta)?,
            rows: CompactRows::default(),
            variable_mod_masses,
            static_mods: view::StaticMods::new(&parameters.static_mods),
            parameters,
            no_protein: Arc::from(""),
            stored_fingerprint: None,
        })
    }

    /// Distinct variable-mod masses, ascending; rows refer to them by index.
    pub fn variable_mod_masses(&self) -> &[f32] {
        &self.variable_mod_masses
    }

    /// The table [`Parameters::digest`] would give, in its row order: generated
    /// directly from the proteome for a non-specific digest
    /// ([`CompactPeptides::build_nonspecific`]), else one length at a time
    /// ([`CompactPeptides::build_by_length`]).
    pub fn build(parameters: &Parameters, fasta: &Fasta) -> Result<Self, CompactError> {
        let enzyme: EnzymeParameters = parameters.enzyme.clone().into();
        match enzyme.enzyme {
            None => Self::build_nonspecific(parameters, fasta),
            Some(_) => Self::build_by_length(parameters, fasta),
        }
    }

    /// Digest one peptide length at a time with sage's own [`Parameters::digest`]
    /// (lengths never interact, so the union is exactly the full digest), compress
    /// each length, then merge the lengths into sage's global row order.
    pub fn build_by_length(parameters: &Parameters, fasta: &Fasta) -> Result<Self, CompactError> {
        let mut table = Self::new(parameters, fasta)?;
        let mut runs = Vec::new();
        for length in table.lengths() {
            let peptides = table.digest_length(fasta, length);
            let run = table.compress(&peptides)?;
            drop_in_parallel(peptides);
            log::info!("compacted {} peptides of length {}", run.len(), length);
            runs.push(run);
        }
        table.merge_length_runs(runs)?;
        Ok(table)
    }

    /// Every peptide length the enzyme settings allow.
    pub fn lengths(&self) -> std::ops::RangeInclusive<usize> {
        let enzyme: EnzymeParameters = self.parameters.enzyme.clone().into();
        enzyme.min_len..=enzyme.max_len
    }

    /// Sage's own digest restricted to peptides of exactly `length` residues.
    pub fn digest_length(&self, fasta: &Fasta, length: usize) -> Vec<Peptide> {
        self.parameters.for_single_length(length).digest(fasta)
    }

    /// Compress `peptides`, keeping their order.
    pub fn compress(&self, peptides: &[Peptide]) -> Result<Vec<CompactPeptide>, CompactError> {
        let target_sequences: Vec<Arc<[u8]>> = peptides
            .par_iter()
            .map(|peptide| match peptide.decoy {
                true => peptide.reverse().sequence,
                false => peptide.sequence.clone(),
            })
            .collect();
        let starts = self.proteome.first_window_starts(&target_sequences);
        drop_in_parallel(target_sequences);

        peptides
            .par_iter()
            .zip(starts)
            .map(|(peptide, start)| self.compress_peptide(peptide, start))
            .collect()
    }

    /// Record as variable-mod sites exactly the slots where `peptide` differs from its
    /// static-mods-only form; [`CompactPeptides::materialize`] replays them.
    fn compress_peptide(&self, peptide: &Peptide, start: u32) -> Result<CompactPeptide, CompactError> {
        let target = match peptide.decoy {
            true => Cow::Owned(peptide.reverse()),
            false => Cow::Borrowed(peptide),
        };
        let mut row = CompactPeptide::new(peptide.monoisotopic, start, &target)?;
        let static_only = self.materialize(&row);
        let unrepresentable = || CompactError::UnrepresentableModification(peptide.to_string());

        let mut applied: Vec<(Site, f32)> = Vec::new();
        for (index, (&mass, &static_mass)) in target
            .modifications
            .iter()
            .zip(&static_only.modifications)
            .enumerate()
        {
            if mass != static_mass {
                applied.push((Site::Sequence(index as u32), mass));
            }
        }
        if target.nterm != static_only.nterm {
            applied.push((Site::Nterm, target.nterm.ok_or_else(unrepresentable)?));
        }
        if target.cterm != static_only.cterm {
            applied.push((Site::Cterm, target.cterm.ok_or_else(unrepresentable)?));
        }
        if applied.len() > MAX_VARIABLE_MODS {
            return Err(CompactError::TooManyVariableMods(applied.len()));
        }
        for (slot, (site, mass)) in row.var_sites.iter_mut().zip(applied) {
            let mass_id = self
                .variable_mod_masses
                .iter()
                .position(|known| known.to_bits() == mass.to_bits())
                .ok_or_else(unrepresentable)?;
            *slot = VarSite::new(site, mass_id as u8);
        }
        if peptide.decoy {
            row.flags |= CompactPeptide::DECOY;
        }
        Ok(row)
    }

    /// The [`Peptide`] `Parameters::digest` holds for `row`, built by the same code
    /// path (`Peptide::try_from(Digest)`, variable sites, static mods, `reverse`).
    /// `proteins` is left empty; see [`CompactPeptides::protein_lists`].
    pub fn materialize(&self, row: &CompactPeptide) -> Peptide {
        let digest = Digest {
            decoy: false,
            semi_enzymatic: row.is_semi_enzymatic(),
            sequence: String::from_utf8(self.proteome.window(row.start, row.length).to_vec())
                .expect("peptide windows are ASCII"),
            protein: self.no_protein.clone(),
            missed_cleavages: row.missed_cleavages(),
            position: row.position(),
        };
        let mut peptide =
            Peptide::try_from(digest).expect("rows only hold sequences sage already accepted");
        peptide.proteins.clear();
        for var_site in row.var_sites {
            if let Some(site) = var_site.site() {
                peptide.apply_site(site, self.variable_mod_masses[var_site.mass_id as usize]);
            }
        }
        for (target, mass) in &self.parameters.static_mods {
            peptide.static_mods(*target, *mass);
        }
        peptide.monoisotopic += peptide.modification_mass();
        match row.is_decoy() {
            true => peptide.reverse(),
            false => peptide,
        }
    }

    /// `materialize` every row in `rows`, in order.
    pub fn materialize_rows(&self, rows: std::ops::Range<usize>) -> Vec<Peptide> {
        rows.into_par_iter()
            .map(|row| self.materialize(&self.rows.get(row)))
            .collect()
    }

    /// How many fragments `build_from_peptides` indexes for `row`: one per ion kind
    /// per backbone cleavage, minus the first `min_ion_index` of each kind.
    pub fn fragment_count(&self, row: &CompactPeptide) -> usize {
        self.fragments_of_length(row.len())
    }

    /// [`CompactPeptides::fragment_count`] of table row `row`.
    pub fn fragment_count_of(&self, row: usize) -> usize {
        self.fragments_of_length(self.rows.length(row))
    }

    fn fragments_of_length(&self, length: usize) -> usize {
        self.parameters.indexed_fragment_count(length)
    }

    /// The peptides a sliced search reports: `rows` (ascending) with their proteins, and
    /// no fragment index. Peptide `i` of it is row `rows[i]`.
    pub fn reported(self: &Arc<Self>, rows: Vec<usize>) -> IndexedDatabase {
        let proteins = self.protein_lists(&rows.iter().map(|&row| self.rows.get(row)).collect::<Vec<_>>());
        IndexedDatabase {
            reported: Some(ReportedPeptides {
                table: Arc::clone(self),
                rows,
                proteins,
            }),
            ion_kinds: self.parameters.ion_kinds.clone(),
            potential_mods: self.parameters.variable_mod_list(),
            bucket_size: self.parameters.bucket_size,
            generate_decoys: self.parameters.generate_decoys,
            decoy_tag: self.parameters.decoy_tag.clone(),
            ..Default::default()
        }
    }

    /// Each row's proteins, sorted, exactly as `Parameters::digest` finds them: every
    /// protein whose own digest yields the row's target form with the same modifications.
    pub fn protein_lists(&self, rows: &[CompactPeptide]) -> Vec<Vec<Arc<str>>> {
        let mut lists = vec![Vec::new(); rows.len()];
        if rows.is_empty() {
            return lists;
        }
        let targets: Vec<CompactView> = rows.iter().map(|&row| self.view_of(row).target()).collect();
        let mut rows_by_sequence: FnvHashMap<&[u8], Vec<usize>> = FnvHashMap::default();
        for (index, target) in targets.iter().enumerate() {
            rows_by_sequence
                .entry(target.target_residues())
                .or_default()
                .push(index);
        }

        let mut enzyme = self.parameters.enzyme.clone();
        enzyme.min_len = targets.iter().map(|t| t.len()).min();
        enzyme.max_len = targets.iter().map(|t| t.len()).max();
        let enzyme: EnzymeParameters = enzyme.into();
        let variable_mods = self.parameters.variable_mod_list();

        let matches: Vec<(usize, usize)> = (0..self.proteome.protein_count())
            .into_par_iter()
            .flat_map_iter(|protein| {
                let accession = self.proteome.accessions[protein].clone();
                let mut found = Vec::new();
                for digest in enzyme.digest(self.proteome.protein_sequence(protein), accession) {
                    let Some(indices) = rows_by_sequence.get(digest.sequence.as_bytes()) else {
                        continue;
                    };
                    let Ok(unmodified) = Peptide::try_from(digest) else {
                        continue;
                    };
                    let variants = unmodified.apply(
                        &variable_mods,
                        &self.parameters.static_mods,
                        self.parameters.max_variable_mods,
                    );
                    for &index in indices {
                        if variants.iter().any(|v| has_same_modified_form(v, &targets[index])) {
                            found.push((index, protein));
                        }
                    }
                }
                found
            })
            .collect();

        for (index, protein) in matches {
            lists[index].push(self.proteome.accessions[protein].clone());
        }
        lists.par_iter_mut().for_each(|list| list.sort_unstable());
        lists
    }

    /// Replace `rows` with the merge of per-length runs (one length per run, each in
    /// its digest's order), giving sage's global row order. Mass ranges merge in parallel.
    pub fn merge_length_runs(&mut self, runs: Vec<Vec<CompactPeptide>>) -> Result<(), CompactError> {
        let total: usize = runs.iter().map(Vec::len).sum();
        if total >= u32::MAX as usize {
            return Err(CompactError::TooManyPeptides(total));
        }
        let pivots = mass_pivots(&runs, rayon::current_num_threads() * 4);
        let bounds: Vec<Vec<usize>> = runs.iter().map(|run| chunk_bounds(run, &pivots)).collect();

        let mut merged = Vec::with_capacity(total);
        merged.par_extend(rayon::iter::repeat_n(CompactPeptide::default(), total));
        let mut outputs = Vec::with_capacity(pivots.len() + 1);
        let mut rest = merged.as_mut_slice();
        for chunk in 0..=pivots.len() {
            let size = bounds.iter().map(|b| b[chunk + 1] - b[chunk]).sum();
            let (output, tail) = std::mem::take(&mut rest).split_at_mut(size);
            outputs.push(output);
            rest = tail;
        }
        outputs
            .into_par_iter()
            .enumerate()
            .for_each(|(chunk, output)| {
                let heads = runs
                    .iter()
                    .zip(&bounds)
                    .map(|(run, b)| &run[b[chunk]..b[chunk + 1]])
                    .collect();
                self.merge_into(heads, output);
            });
        self.rows = CompactRows::from_rows(merged);
        Ok(())
    }

    /// `a` (a row of this table) and `b` (a row of `other`, built over the same
    /// FASTA) are the same peptide, read from the same or another occurrence of its
    /// residues: equal length, decoy flag, variable sites and residues.
    pub fn same_row_apart_from_occurrence(
        &self,
        a: &CompactPeptide,
        other: &CompactPeptides,
        b: &CompactPeptide,
    ) -> bool {
        a.length == b.length
            && a.is_decoy() == b.is_decoy()
            && a.var_sites == b.var_sites
            && self.proteome.window(a.start, a.length) == other.proteome.window(b.start, b.length)
    }

    /// Adjacent rows out of sage's order. Rows of different lengths are compared in
    /// full; rows of one length keep their digest's order, so only by mass.
    pub fn row_order_violations(&self) -> usize {
        (1..self.rows.len())
            .into_par_iter()
            .filter(|&next| {
                let (a, b) = (&self.rows.get(next - 1), &self.rows.get(next));
                match a.length == b.length {
                    true => a.monoisotopic.total_cmp(&b.monoisotopic).is_gt(),
                    false => !self.compare_different_lengths(a, b).is_lt(),
                }
            })
            .count()
    }

    /// K-way merge of runs that each hold a single, distinct peptide length.
    fn merge_into(&self, mut heads: Vec<&[CompactPeptide]>, output: &mut [CompactPeptide]) {
        for slot in output.iter_mut() {
            let mut best: Option<usize> = None;
            for (run, head) in heads.iter().enumerate() {
                let Some(candidate) = head.first() else {
                    continue;
                };
                let better = match best {
                    None => true,
                    Some(b) => self.compare_different_lengths(candidate, &heads[b][0]).is_lt(),
                };
                if better {
                    best = Some(run);
                }
            }
            let best = best.expect("output is sized to the rows the heads hold");
            *slot = heads[best][0];
            heads[best] = &heads[best][1..];
        }
    }

    /// Sage's row order (`Parameters::reorder_peptides`: mass, then
    /// `Peptide::initial_sort`) for two peptides of different lengths, whose
    /// sequences can never be equal, so modifications are never consulted.
    fn compare_different_lengths(&self, a: &CompactPeptide, b: &CompactPeptide) -> Ordering {
        a.monoisotopic
            .total_cmp(&b.monoisotopic)
            .then_with(|| self.compare_sequences(a, b))
    }

    fn compare_sequences(&self, a: &CompactPeptide, b: &CompactPeptide) -> Ordering {
        let shared = a.len().min(b.len());
        for index in 0..shared {
            let order = self.residue(a, index).cmp(&self.residue(b, index));
            if order.is_ne() {
                return order;
            }
        }
        a.len().cmp(&b.len())
    }

    /// Residue `index` of `row`'s sequence; a decoy reads its target's interior
    /// backwards, as [`Peptide::reverse`] builds it.
    fn residue(&self, row: &CompactPeptide, index: usize) -> u8 {
        let length = row.len();
        let in_reversed_interior = row.is_decoy() && length > 2 && index > 0 && index < length - 1;
        let target_index = match in_reversed_interior {
            true => length - 1 - index,
            false => index,
        };
        self.proteome.residues[row.start as usize + target_index]
    }
}

/// Free a large vector's elements from all threads. Dropping tens of millions of
/// small heap allocations on one thread took longer than building them, because
/// glibc serialises frees of memory other threads allocated.
pub fn drop_in_parallel<T: Send>(items: Vec<T>) {
    items.into_par_iter().for_each(drop);
}

fn has_same_modified_form(a: &impl PeptideView, b: &impl PeptideView) -> bool {
    a.len() == b.len()
        && (0..a.len()).all(|i| a.residue(i) == b.residue(i) && a.modification(i) == b.modification(i))
        && a.nterm() == b.nterm()
        && a.cterm() == b.cterm()
}

/// Ascending, distinct masses that split all runs together into about `chunks`
/// equally sized parts, from a sample of every run (runs of different lengths
/// cover different mass ranges).
fn mass_pivots(runs: &[Vec<CompactPeptide>], chunks: usize) -> Vec<f32> {
    const SAMPLE_EVERY: usize = 1024;
    let mut sample: Vec<f32> = runs
        .iter()
        .flat_map(|run| run.iter().step_by(SAMPLE_EVERY).map(|row| row.monoisotopic))
        .collect();
    if sample.is_empty() {
        return Vec::new();
    }
    sample.sort_unstable_by(f32::total_cmp);
    let mut pivots: Vec<f32> = (1..chunks)
        .map(|chunk| sample[chunk * sample.len() / chunks])
        .collect();
    pivots.dedup_by(|a, b| a.to_bits() == b.to_bits());
    pivots
}

/// `[0, rows below pivot 0, rows below pivot 1, ..., run.len()]`.
fn chunk_bounds(run: &[CompactPeptide], pivots: &[f32]) -> Vec<usize> {
    let mut bounds = Vec::with_capacity(pivots.len() + 2);
    bounds.push(0);
    for pivot in pivots {
        bounds.push(run.partition_point(|row| row.monoisotopic.total_cmp(pivot).is_lt()));
    }
    bounds.push(run.len());
    bounds
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::database::{Builder, PeptideIx};

    const FASTA: &str = "\
>sp|P1|ONE
MQPEPTIDEKCMAMRPLLGQKSTWCMKAAMMK
>sp|P2|TWO
PEPTIDEKCMAMRAAQMGHKCPEPTIDEK
>sp|P3|THREE
QMGHKLLWPEPTIDEKMMCQ
>rev_sp|P4|DECOY
AAAAKCMR
>sp|P5|FIVE
MCQPEPTIDEK
";

    type DuplicateAttributes = (Position, bool, u8);

    fn duplicate_attributes(peptide: &Peptide) -> DuplicateAttributes {
        (peptide.position, peptide.semi_enzymatic, peptide.missed_cleavages)
    }

    /// A peptide merged from several digests keeps the attributes of whichever
    /// duplicate `reorder_peptides`' unstable sort happened to keep, which depends on
    /// everything else in the array.
    fn without_duplicate_attributes(peptide: &Peptide) -> Peptide {
        Peptide {
            position: Position::Internal,
            semi_enzymatic: false,
            missed_cleavages: 0,
            ..peptide.clone()
        }
    }

    /// The attributes of every digest that yields `peptide`'s modified target form.
    fn attributes_of_duplicates(
        parameters: &Parameters,
        fasta: &Fasta,
        peptide: &Peptide,
    ) -> Vec<DuplicateAttributes> {
        let target = match peptide.decoy {
            true => peptide.reverse(),
            false => peptide.clone(),
        };
        let enzyme: EnzymeParameters = parameters.enzyme.clone().into();
        let mut attributes = Vec::new();
        for (accession, sequence) in &fasta.targets {
            for digest in enzyme.digest(sequence, accession.clone()) {
                if digest.sequence.as_bytes() != &target.sequence[..] {
                    continue;
                }
                let variants = Peptide::try_from(digest).unwrap().apply(
                    &parameters.variable_mod_list(),
                    &parameters.static_mods,
                    parameters.max_variable_mods,
                );
                for variant in variants {
                    if has_same_modified_form(&variant, &target) {
                        attributes.push(duplicate_attributes(&variant));
                    }
                }
            }
        }
        attributes
    }

    fn parameters(config: serde_json::Value) -> Parameters {
        let mut builder: Builder = serde_json::from_value(config).unwrap();
        builder.fasta = Some("unused".into());
        builder.make_parameters()
    }

    fn assert_matches_sage_digest(config: serde_json::Value) {
        assert_fasta_matches_sage_digest(FASTA, config);
    }

    fn assert_fasta_matches_sage_digest(fasta: &str, config: serde_json::Value) {
        let fasta = Fasta::parse(fasta.into(), "rev_", true);
        let parameters = parameters(config);
        let expected = parameters.digest(&fasta);
        assert!(!expected.is_empty());

        let table = CompactPeptides::build(&parameters, &fasta).unwrap();
        let mut actual: Vec<Peptide> = table
            .rows
            .iter()
            .map(|row| table.materialize(&row))
            .collect();
        let rows: Vec<CompactPeptide> = table.rows.iter().collect();
        for (peptide, proteins) in actual.iter_mut().zip(table.protein_lists(&rows)) {
            peptide.proteins = proteins;
        }

        assert_eq!(actual.len(), expected.len());
        for (row, ((actual, expected), compact)) in actual
            .iter()
            .zip(&expected)
            .zip(table.rows.iter())
            .enumerate()
        {
            assert_eq!(
                compact.monoisotopic.to_bits(),
                expected.monoisotopic.to_bits(),
                "row {row}: compact mass {} vs sage {expected}",
                compact.monoisotopic
            );
            assert_eq!(
                actual.monoisotopic.to_bits(),
                expected.monoisotopic.to_bits(),
                "row {row}: materialized mass vs sage {expected}"
            );
            assert_eq!(
                without_duplicate_attributes(actual),
                without_duplicate_attributes(expected),
                "row {row}: compact {actual} vs sage {expected}"
            );
            let allowed = attributes_of_duplicates(&parameters, &fasta, expected);
            assert!(
                allowed.contains(&duplicate_attributes(expected)),
                "row {row}: sage's own {expected} attributes not among its duplicates {allowed:?}"
            );
            assert!(
                allowed.contains(&duplicate_attributes(actual)),
                "row {row}: compact {actual} attributes {:?} not among its duplicates {allowed:?}",
                duplicate_attributes(actual)
            );
            assert_eq!(
                actual.modifications, expected.modifications,
                "row {row}: {actual} vs {expected}"
            );
        }
    }

    fn mods(config: serde_json::Value) -> serde_json::Value {
        let mut base = serde_json::json!({
            "static_mods": {"C": 57.021465},
            "variable_mods": {"M": [15.9949], "[": [42.010565], "^Q": [-17.026548]},
            "max_variable_mods": 2,
            "peptide_min_mass": 0.0,
            "peptide_max_mass": 5000.0,
            "generate_decoys": true,
            "decoy_tag": "rev_",
        });
        base.as_object_mut()
            .unwrap()
            .extend(config.as_object().unwrap().clone());
        base
    }

    #[test]
    fn tryptic_matches_sage_digest() {
        assert_matches_sage_digest(mods(serde_json::json!({
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 2, "min_len": 3, "max_len": 30},
        })));
    }

    #[test]
    fn semi_tryptic_matches_sage_digest() {
        assert_matches_sage_digest(mods(serde_json::json!({
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 1, "min_len": 4, "max_len": 20, "semi_enzymatic": true},
        })));
    }

    #[test]
    fn nonspecific_matches_sage_digest() {
        assert_matches_sage_digest(mods(serde_json::json!({
            "enzyme": {"cleave_at": "", "min_len": 4, "max_len": 9},
        })));
    }

    #[test]
    fn nonspecific_terminal_mods_match_sage_digest() {
        assert_matches_sage_digest(mods(serde_json::json!({
            "enzyme": {"cleave_at": "", "min_len": 5, "max_len": 8},
            "static_mods": {"^": 229.162932, "C": 57.021465},
            "variable_mods": {"M": [15.9949, 31.989829], "[": [42.010565], "$": [-0.984016], "^Q": [-17.026548]},
        })));
    }

    #[test]
    fn immunopeptidomics_matches_sage_digest() {
        assert_matches_sage_digest(serde_json::json!({
            "enzyme": {"cleave_at": "", "min_len": 6, "max_len": 10},
            "static_mods": {},
            "variable_mods": {"M": [15.9949], "[": [42.010565], "C": [119.004099]},
            "max_variable_mods": 1,
            "peptide_min_mass": 0.0,
            "peptide_max_mass": 5000.0,
            "generate_decoys": true,
            "decoy_tag": "rev_",
        }));
    }

    /// Non-specific edge cases: a C-terminal window repeating an internal one, an
    /// N-terminal window repeated later (also as the C-terminal one), an `X`, a
    /// protein shorter than any window, a decoy equal to a target window (PEMKQ
    /// reversed inside is PKMEQ), a palindromic interior, a decoy-tagged protein, a
    /// protein duplicated whole and a whole-protein window.
    const EDGE_FASTA: &str = "\
>sp|E1|CTERM_REPEAT
AGKLMNPWQAKLMNP
>sp|E2|NTERM_REPEAT
MSTAYGGHKMSTAY
>sp|E3|WITH_X
PEPXTIDEKAMCQ
>sp|E4|SHORT
MK
>sp|E5|DECOY_EQUALS_TARGET
PEMKQRSTPKMEQ
>sp|E6|PALINDROME
GAKMKAG
>rev_sp|E7|DECOY
AAAAKCMR
>sp|E8|WHOLE_DUPLICATE
MSTAYGGHKMSTAY
>sp|E9|FULL
MQKCPE
";

    fn nonspecific_configs() -> Vec<serde_json::Value> {
        vec![
            mods(serde_json::json!({
                "enzyme": {"cleave_at": "", "min_len": 4, "max_len": 8},
            })),
            mods(serde_json::json!({
                "enzyme": {"cleave_at": "", "min_len": 4, "max_len": 8},
                "static_mods": {"[": 42.010565, "C": 57.021465},
                "variable_mods": {"M": [15.9949], "]": [-0.984016], "^Q": [-17.026548]},
            })),
            mods(serde_json::json!({
                "enzyme": {"cleave_at": "", "min_len": 3, "max_len": 9},
                "peptide_min_mass": 500.0,
                "peptide_max_mass": 800.0,
            })),
            serde_json::json!({
                "enzyme": {"cleave_at": "", "min_len": 4, "max_len": 9},
                "static_mods": {},
                "variable_mods": {"M": [15.9949], "[": [42.010565], "C": [119.004099]},
                "max_variable_mods": 1,
                "peptide_min_mass": 0.0,
                "peptide_max_mass": 5000.0,
                "generate_decoys": true,
                "decoy_tag": "rev_",
            }),
        ]
    }

    #[test]
    fn nonspecific_edge_cases_match_sage_digest() {
        for config in nonspecific_configs() {
            assert_fasta_matches_sage_digest(EDGE_FASTA, config);
        }
    }

    /// The same rows in the same order; `start` may point at another occurrence of
    /// the same residues.
    #[test]
    fn nonspecific_build_equals_build_by_length() {
        for fasta in [FASTA, EDGE_FASTA] {
            for config in nonspecific_configs() {
                let fasta = Fasta::parse(fasta.into(), "rev_", true);
                let parameters = parameters(config);
                let direct = CompactPeptides::build_nonspecific(&parameters, &fasta).unwrap();
                let by_length = CompactPeptides::build_by_length(&parameters, &fasta).unwrap();
                assert_eq!(direct.rows.len(), by_length.rows.len());
                for (row, (a, b)) in direct.rows.iter().zip(by_length.rows.iter()).enumerate() {
                    assert_eq!(a.monoisotopic.to_bits(), b.monoisotopic.to_bits(), "row {row}");
                    assert_eq!(
                        (a.length, a.is_decoy(), a.var_sites),
                        (b.length, b.is_decoy(), b.var_sites),
                        "row {row}"
                    );
                    assert_eq!(
                        direct.proteome.window(a.start, a.length),
                        by_length.proteome.window(b.start, b.length),
                        "row {row}"
                    );
                }
            }
        }
    }

    fn tryptic() -> serde_json::Value {
        mods(serde_json::json!({
            "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 2, "min_len": 3, "max_len": 30},
        }))
    }

    #[test]
    fn fragment_count_matches_the_built_index() {
        for config in [tryptic(), mods(serde_json::json!({"enzyme": {"cleave_at": "", "min_len": 4, "max_len": 9}}))] {
            let fasta = Fasta::parse(FASTA.into(), "rev_", true);
            let parameters = parameters(config);
            let table = CompactPeptides::build(&parameters, &fasta).unwrap();
            let counted: usize = table
                .rows
                .iter()
                .map(|row| table.fragment_count(&row))
                .sum();
            let built = parameters
                .clone()
                .build_from_peptides(parameters.digest(&fasta))
                .size();
            assert_eq!(counted, built);
        }
    }

    #[test]
    fn reported_peptides_read_as_sage_digests_them() {
        let fasta = Fasta::parse(FASTA.into(), "rev_", true);
        let parameters = parameters(tryptic());
        let expected = parameters.digest(&fasta);
        let table = Arc::new(CompactPeptides::build(&parameters, &fasta).unwrap());
        let rows: Vec<usize> = (0..table.rows.len()).step_by(3).collect();
        let database = table.reported(rows.clone());
        assert_eq!(database.row_count(), rows.len());
        for (index, &row) in rows.iter().enumerate() {
            let (peptide, expected) = (database.peptide(PeptideIx(index as u32)), &expected[row]);
            assert_eq!(peptide.to_string(), expected.to_string());
            assert_eq!(peptide.target_string(), match expected.decoy {
                true => expected.reverse().to_string(),
                false => expected.to_string(),
            });
            assert_eq!(peptide.decoy(), expected.decoy);
            assert_eq!(database.proteins(PeptideIx(index as u32)), &expected.proteins[..]);
            assert_eq!(
                database.protein_names(PeptideIx(index as u32)),
                expected.proteins("rev_", true)
            );
        }
        assert_eq!(database.decoy_tag, "rev_");
        assert!(database.fragment_mzs.is_empty() && database.peptides.is_empty());
    }

    #[test]
    fn rejects_more_than_two_variable_mods() {
        let fasta = Fasta::parse(FASTA.into(), "rev_", true);
        let parameters = parameters(mods(serde_json::json!({"max_variable_mods": 3})));
        assert_eq!(
            CompactPeptides::new(&parameters, &fasta).err(),
            Some(CompactError::TooManyVariableMods(3))
        );
    }

    #[test]
    fn rows_are_sixteen_bytes() {
        assert_eq!(std::mem::size_of::<CompactPeptide>(), 16);
    }

    /// Every row's view reads as its materialised `Peptide` does: residues,
    /// modification bits, termini, mass, decoy, missed cleavages and every ion.
    fn assert_views_equal_materialized(fasta: &str, config: serde_json::Value) {
        let fasta = Fasta::parse(fasta.into(), "rev_", true);
        let parameters = parameters(config.clone());
        let table = CompactPeptides::build(&parameters, &fasta).unwrap();
        assert!(!table.rows.is_empty());
        for row in 0..table.rows.len() {
            let peptide = table.materialize(&table.rows.get(row));
            if let Some(difference) =
                crate::peptide_view::first_difference(&table.view(row), &peptide)
            {
                panic!("row {row} ({peptide}): {difference}; config {config}");
            }
        }
    }

    /// Static rules that compete for the same residue or terminus, with variable
    /// sites on the same residues: which one sticks depends on the order they apply.
    fn colliding_mods() -> serde_json::Value {
        mods(serde_json::json!({
            "enzyme": {"cleave_at": "", "min_len": 3, "max_len": 8},
            "static_mods": {
                "C": 57.021465, "^C": 10.5, "^M": 3.25, "[M": 42.010565, "M": 1.75,
                "$K": 1.5, "]K": 2.5, "K": 0.125, "^": 229.162932, "]": -0.984016
            },
            "variable_mods": {"M": [15.9949], "C": [119.004099], "[": [42.010565], "^Q": [-17.026548], "$": [14.01565]},
        }))
    }

    /// Protein-terminal residue rules with nothing competing, so only the row's
    /// position decides whether they apply.
    fn protein_terminal_mods() -> serde_json::Value {
        mods(serde_json::json!({
            "enzyme": {"cleave_at": "", "min_len": 3, "max_len": 8},
            "static_mods": {"[P": 7.5, "[A": 6.5, "]Q": 4.25, "]P": 3.75},
        }))
    }

    /// An index over compact rows holds what the same rows materialised give it: row
    /// masses, fragments and every query's precursor range.
    #[test]
    fn index_over_compact_rows_equals_index_over_materialized_rows() {
        use crate::fragment_grid::{FragmentLayout, GridShape};
        use crate::mass::Tolerance;
        let fasta = Fasta::parse(FASTA.into(), "rev_", true);
        let grid = GridShape {
            rows_per_block: 7,
            fine_bits: 10,
            cell_fragments: 2,
        };
        for layout in [FragmentLayout::Pages, FragmentLayout::Grid(grid)] {
            let mut parameters = parameters(colliding_mods());
            parameters.fragment_layout = layout;
            parameters.bucket_size = 64;
            let table = Arc::new(CompactPeptides::build(&parameters, &fasta).unwrap());
            let rows = 3..table.rows.len() - 2;
            let peptides = table.materialize_rows(rows.clone());
            let built = parameters.clone().build_from_peptides_with_rt(peptides, None);
            let compact = parameters
                .clone()
                .build_from_compact(Arc::clone(&table), rows.clone(), None);
            assert!(compact.peptides.is_empty());
            assert_eq!(built.row_count(), compact.row_count());
            for row in 0..built.row_count() {
                assert_eq!(built.row_mass(row).to_bits(), compact.row_mass(row).to_bits());
            }
            let fragments = |db: &IndexedDatabase| {
                let mut fragments: Vec<(u32, u32)> = db
                    .fragments()
                    .map(|f| (f.peptide_index.0, f.fragment_mz.to_bits()))
                    .collect();
                fragments.sort_unstable();
                fragments
            };
            assert_eq!(fragments(&built), fragments(&compact));
            for row in 0..built.row_count() {
                let mass = built.row_mass(row);
                for tolerance in [Tolerance::Ppm(-10.0, 10.0), Tolerance::Da(-1.5, 0.5)] {
                    let (a, b) = (built.query(mass, tolerance), compact.query(mass, tolerance));
                    assert_eq!((a.pre_idx_lo, a.pre_idx_hi), (b.pre_idx_lo, b.pre_idx_hi));
                    assert_eq!(a.candidate_rows(), b.candidate_rows());
                }
            }
        }
    }

    #[test]
    fn views_equal_materialized_peptides() {
        for config in [
            mods(serde_json::json!({
                "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 2, "min_len": 3, "max_len": 30},
            })),
            mods(serde_json::json!({
                "enzyme": {"cleave_at": "KR", "restrict": "P", "missed_cleavages": 1, "min_len": 4, "max_len": 20, "semi_enzymatic": true},
            })),
            mods(serde_json::json!({
                "enzyme": {"cleave_at": "", "min_len": 5, "max_len": 8},
                "static_mods": {"^": 229.162932, "C": 57.021465},
                "variable_mods": {"M": [15.9949, 31.989829], "[": [42.010565], "$": [-0.984016], "^Q": [-17.026548]},
            })),
            colliding_mods(),
            protein_terminal_mods(),
        ] {
            assert_views_equal_materialized(FASTA, config);
        }
        for config in nonspecific_configs()
            .into_iter()
            .chain([colliding_mods(), protein_terminal_mods()])
        {
            assert_views_equal_materialized(EDGE_FASTA, config);
        }
    }
}
