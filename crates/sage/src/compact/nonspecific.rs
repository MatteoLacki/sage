//! [`CompactPeptides`] for a non-specific digest, generated straight from the
//! proteome without building any `Peptide`. Design and verification:
//! `docs/ai/compact_peptides.md`.

use super::{
    CompactError, CompactPeptide, CompactPeptides, Proteome, VarSite, MAX_PEPTIDE_LENGTH,
    MAX_VARIABLE_MODS,
};
use crate::database::Parameters;
use crate::enzyme::{EnzymeParameters, Position};
use crate::fasta::Fasta;
use crate::mass::{monoisotopic, H2O};
use crate::modification::ModificationSpecificity;
use crate::peptide::Site;
use fnv::FnvHasher;
use rayon::prelude::*;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hasher;
use std::marker::PhantomData;
use std::ops::{Range, RangeInclusive};
use std::sync::atomic::{AtomicU32, AtomicU64, AtomicUsize, Ordering as AtomicOrdering};
use std::time::Instant;

/// Width of the mass bins rows are counted by.
const FINE_BIN_DA: f64 = 0.01;
const MAX_FINE_BINS: usize = 1 << 22;
/// Rows per placement bin; each bin is sorted on its own, in about 1 MB.
const ROWS_PER_PLACEMENT_BIN: u64 = 1 << 16;
/// Protein chunks per thread, so chunks of uneven cost still balance.
const CHUNKS_PER_THREAD: usize = 8;

impl CompactPeptides {
    /// [`Parameters::digest`]'s table for a non-specific digest, generated straight
    /// from the proteome: rows are counted per mass bin, generated again straight
    /// into their bins, then each bin is sorted into sage's order and its duplicates
    /// merged.
    pub fn build_nonspecific(parameters: &Parameters, fasta: &Fasta) -> Result<Self, CompactError> {
        let mut table = Self::new(parameters, fasta)?;
        let digest = NonspecificDigest::new(parameters);

        let phase = Instant::now();
        let decoys = DecoyFilter::build(&table.proteome, digest.lengths());
        let filter_seconds = phase.elapsed().as_secs_f64();

        let phase = Instant::now();
        let chunks = protein_chunks(
            &table.proteome,
            rayon::current_num_threads() * CHUNKS_PER_THREAD,
        );
        let mass_bins = MassBins::new(digest.min_mass, digest.max_mass);
        let counts: Vec<AtomicU64> = (0..mass_bins.count)
            .into_par_iter()
            .map(|_| AtomicU64::new(0))
            .collect();
        chunks.par_iter().try_for_each(|proteins| {
            table.generate_rows(&digest, &decoys, proteins.clone(), &mut |row| {
                counts[mass_bins.of(row.monoisotopic)].fetch_add(1, AtomicOrdering::Relaxed);
            })
        })?;
        let counts: Vec<u64> = counts.into_iter().map(AtomicU64::into_inner).collect();
        let placement = Placement::new(&counts);
        let count_seconds = phase.elapsed().as_secs_f64();

        let phase = Instant::now();
        let mut rows = Vec::with_capacity(placement.total());
        rows.par_extend(rayon::iter::repeat_n(
            CompactPeptide::default(),
            placement.total(),
        ));
        let cursors: Vec<AtomicUsize> = (0..placement.bins())
            .map(|bin| AtomicUsize::new(placement.range(bin).start))
            .collect();
        let slots = RowSlots::new(&mut rows);
        chunks.par_iter().try_for_each(|proteins| {
            table.generate_rows(&digest, &decoys, proteins.clone(), &mut |row| {
                let bin = placement.bin_of_fine[mass_bins.of(row.monoisotopic)] as usize;
                let index = cursors[bin].fetch_add(1, AtomicOrdering::Relaxed);
                assert!(
                    index < placement.range(bin).end,
                    "placing generated more rows into a mass bin than counting did"
                );
                // SAFETY: `fetch_add` hands each index of the bin to exactly one caller,
                // and the assert keeps it inside the bin.
                unsafe { slots.write(index, row) };
            })
        })?;
        for (bin, cursor) in cursors.iter().enumerate() {
            assert_eq!(
                cursor.load(AtomicOrdering::Relaxed),
                placement.range(bin).end,
                "placing generated fewer rows into a mass bin than counting did"
            );
        }
        let place_seconds = phase.elapsed().as_secs_f64();

        let phase = Instant::now();
        let survivors = table.sort_and_merge_bins(&mut rows, &placement);
        let mut packed = 0;
        for (bin, &kept) in survivors.iter().enumerate() {
            let from = placement.range(bin).start;
            rows.copy_within(from..from + kept, packed);
            packed += kept;
        }
        rows.truncate(packed);
        rows.shrink_to_fit();
        if rows.len() >= u32::MAX as usize {
            return Err(CompactError::TooManyPeptides(rows.len()));
        }
        log::info!(
            "non-specific digest: decoy filter {:.1} s, {} rows counted in {:.1} s, placed in {:.1} s, \
             sorted and merged into {} in {:.1} s",
            filter_seconds,
            placement.total(),
            count_seconds,
            place_seconds,
            rows.len(),
            phase.elapsed().as_secs_f64(),
        );
        table.rows = rows;
        Ok(table)
    }

    /// Every row [`Parameters::digest`] yields from `proteins`, targets and kept
    /// decoys, passed to `emit` in no particular order. A sequence occurring several
    /// times is emitted from each occurrence that yields the same rows as its first
    /// one; [`CompactPeptides::sort_and_merge_bins`] merges them.
    fn generate_rows(
        &self,
        digest: &NonspecificDigest,
        decoys: &DecoyFilter,
        proteins: Range<usize>,
        emit: &mut impl FnMut(CompactPeptide),
    ) -> Result<(), CompactError> {
        let mut scratch = Scratch {
            variable_sites: Vec::new(),
            static_only: ModifiedForm::new(),
            form: ModifiedForm::new(),
        };
        for protein in proteins {
            let range = self.proteome.protein_range(protein);
            let residues = &self.proteome.residues[range.clone()];
            for first in 0..residues.len() {
                // `Peptide::try_from`'s mass: `H2O`, then each residue in order.
                let mut mass = H2O;
                let longest = digest.max_len.min(residues.len() - first);
                for length in 1..=longest {
                    let residue_mass = monoisotopic(residues[first + length - 1]);
                    if residue_mass == 0.0 {
                        break;
                    }
                    mass += residue_mass;
                    if length < digest.min_len {
                        continue;
                    }
                    let position = window_position(first, length, residues.len());
                    if repeats_an_occurrence_with_another_position(
                        residues, first, length, position,
                    ) {
                        continue;
                    }
                    let start = (range.start + first) as u32;
                    let keep_decoy = !decoys.is_dropped(start as usize, length);
                    self.generate_window_rows(
                        digest,
                        &residues[first..first + length],
                        start,
                        position,
                        mass,
                        keep_decoy,
                        &mut scratch,
                        emit,
                    )?;
                }
            }
        }
        Ok(())
    }

    /// The rows `Peptide::apply` and the decoy step of [`Parameters::digest`] make of
    /// one window: the unmodified form and every combination of variable sites, each
    /// with static mods, inside the mass range.
    #[allow(clippy::too_many_arguments)]
    fn generate_window_rows(
        &self,
        digest: &NonspecificDigest,
        residues: &[u8],
        start: u32,
        position: Position,
        mass: f32,
        keep_decoy: bool,
        scratch: &mut Scratch,
        emit: &mut impl FnMut(CompactPeptide),
    ) -> Result<(), CompactError> {
        let Scratch {
            variable_sites,
            static_only,
            form,
        } = scratch;
        static_only.reset(residues.len());
        static_only.apply_static_mods(&self.parameters.static_mods, residues, position);
        variable_sites.clear();
        for &(target, mod_mass) in &digest.variable_mods {
            for_each_site(target, residues, position, |site| {
                variable_sites.push((site, mod_mass))
            });
        }

        let mut emit_combination = |combination: &[(Site, f32)]| -> Result<(), CompactError> {
            form.reset(residues.len());
            for &(site, mod_mass) in combination {
                form.apply_site(site, mod_mass);
            }
            form.apply_static_mods(&self.parameters.static_mods, residues, position);
            let monoisotopic = mass + form.modification_mass();
            if !(monoisotopic >= digest.min_mass && monoisotopic <= digest.max_mass) {
                return Ok(());
            }
            let var_sites = self.variable_sites_against(form, static_only, residues)?;
            let row = CompactPeptide::nonspecific_target(
                monoisotopic,
                start,
                residues.len(),
                position,
                var_sites,
            );
            emit(row);
            if keep_decoy {
                emit(row.decoy_of());
            }
            Ok(())
        };

        emit_combination(&[])?;
        if digest.variable_mods.is_empty() {
            return Ok(());
        }
        let sites: &[(Site, f32)] = variable_sites;
        if digest.max_variable_mods >= 1 {
            for site in sites {
                emit_combination(std::slice::from_ref(site))?;
            }
        }
        if digest.max_variable_mods >= 2 {
            for (index, &a) in sites.iter().enumerate() {
                for &b in &sites[index + 1..] {
                    if a.0 != b.0 {
                        emit_combination(&[a, b])?;
                    }
                }
            }
        }
        Ok(())
    }

    /// The variable sites `compress_peptide` records for `form`: the slots where it
    /// differs from `static_only`, sequence sites ascending, then N- and C-terminus.
    fn variable_sites_against(
        &self,
        form: &ModifiedForm,
        static_only: &ModifiedForm,
        residues: &[u8],
    ) -> Result<[VarSite; MAX_VARIABLE_MODS], CompactError> {
        let unrepresentable = || {
            CompactError::UnrepresentableModification(
                String::from_utf8_lossy(residues).into_owned(),
            )
        };
        let mut sites = [VarSite::UNUSED; MAX_VARIABLE_MODS];
        let mut count = 0;
        let mut record = |site: Site, mass: f32| -> Result<(), CompactError> {
            if count < MAX_VARIABLE_MODS {
                let mass_id = self
                    .variable_mod_masses
                    .iter()
                    .position(|known| known.to_bits() == mass.to_bits())
                    .ok_or_else(unrepresentable)?;
                sites[count] = VarSite::new(site, mass_id as u8);
            }
            count += 1;
            Ok(())
        };
        for (index, (&mass, &static_mass)) in form
            .modifications()
            .iter()
            .zip(static_only.modifications())
            .enumerate()
        {
            if mass != static_mass {
                record(Site::Sequence(index as u32), mass)?;
            }
        }
        if form.nterm != static_only.nterm {
            record(Site::Nterm, form.nterm.ok_or_else(unrepresentable)?)?;
        }
        if form.cterm != static_only.cterm {
            record(Site::Cterm, form.cterm.ok_or_else(unrepresentable)?)?;
        }
        if count > MAX_VARIABLE_MODS {
            return Err(CompactError::TooManyVariableMods(count));
        }
        Ok(sites)
    }

    /// Sort each placement bin into sage's order and merge its duplicates in place.
    /// Returns each bin's survivor count; survivors are packed at the bin's front.
    fn sort_and_merge_bins(
        &self,
        rows: &mut [CompactPeptide],
        placement: &Placement,
    ) -> Vec<usize> {
        let mut bins = Vec::with_capacity(placement.bins());
        let mut rest = rows;
        for bin in 0..placement.bins() {
            let (head, tail) = std::mem::take(&mut rest).split_at_mut(placement.range(bin).len());
            bins.push(head);
            rest = tail;
        }
        bins.into_par_iter()
            .map(|bin| {
                // Among equal peptides the smallest start survives: the first
                // occurrence in the first protein holding it.
                bin.sort_unstable_by(|a, b| {
                    self.compare_rows(a, b)
                        .then_with(|| a.start.cmp(&b.start))
                        .then_with(|| a.flags.cmp(&b.flags))
                });
                self.merge_duplicates(bin)
            })
            .collect()
    }

    /// `Parameters::reorder_peptides`' `dedup_by` on a sorted run: the first of each
    /// group of equal peptides survives. Returns the survivor count.
    fn merge_duplicates(&self, rows: &mut [CompactPeptide]) -> usize {
        if rows.is_empty() {
            return 0;
        }
        let mut kept = 1;
        for next in 1..rows.len() {
            if self.same_peptide(&rows[kept - 1], &rows[next]) {
                debug_assert_eq!(rows[kept - 1].is_decoy(), rows[next].is_decoy());
                continue;
            }
            rows[kept] = rows[next];
            kept += 1;
        }
        kept
    }

    /// `Parameters::reorder_peptides`' order: mass, then `Peptide::initial_sort`
    /// (sequence, modifications, N-terminus, C-terminus).
    fn compare_rows(&self, a: &CompactPeptide, b: &CompactPeptide) -> Ordering {
        a.monoisotopic
            .total_cmp(&b.monoisotopic)
            .then_with(|| self.compare_sequences(a, b))
            .then_with(|| self.compare_modified_forms(a, b))
    }

    /// `Parameters::reorder_peptides`' dedup condition.
    fn same_peptide(&self, a: &CompactPeptide, b: &CompactPeptide) -> bool {
        if a.monoisotopic != b.monoisotopic || self.compare_sequences(a, b).is_ne() {
            return false;
        }
        if have_the_same_form_by_construction(a, b) {
            return true;
        }
        let (a, b) = (self.modified_form(a), self.modified_form(b));
        a.modifications() == b.modifications() && a.nterm == b.nterm && a.cterm == b.cterm
    }

    /// `Peptide::initial_sort` past the sequence, for two rows of equal sequence.
    fn compare_modified_forms(&self, a: &CompactPeptide, b: &CompactPeptide) -> Ordering {
        if have_the_same_form_by_construction(a, b) {
            return Ordering::Equal;
        }
        let (a, b) = (self.modified_form(a), self.modified_form(b));
        a.modifications()
            .partial_cmp(b.modifications())
            .unwrap_or(Ordering::Equal)
            .then_with(|| a.nterm.partial_cmp(&b.nterm).unwrap_or(Ordering::Equal))
            .then_with(|| a.cterm.partial_cmp(&b.cterm).unwrap_or(Ordering::Equal))
    }

    /// The `modifications`, `nterm` and `cterm` [`CompactPeptides::materialize`]
    /// gives `row`.
    fn modified_form(&self, row: &CompactPeptide) -> ModifiedForm {
        let residues = self.proteome.window(row.start, row.length);
        let mut form = ModifiedForm::new();
        form.reset(row.len());
        for var_site in row.var_sites {
            if let Some(site) = var_site.site() {
                form.apply_site(site, self.variable_mod_masses[var_site.mass_id as usize]);
            }
        }
        form.apply_static_mods(&self.parameters.static_mods, residues, row.position());
        if row.is_decoy() {
            form.reverse_interior();
        }
        form
    }
}

/// Two rows of equal sequence: a row's modified form is fixed by its target residues,
/// position, variable sites and decoy flag.
fn have_the_same_form_by_construction(a: &CompactPeptide, b: &CompactPeptide) -> bool {
    a.var_sites == b.var_sites && a.position() == b.position() && a.is_decoy() == b.is_decoy()
}

/// The digest settings the generator reads.
struct NonspecificDigest {
    min_len: usize,
    max_len: usize,
    min_mass: f32,
    max_mass: f32,
    variable_mods: Vec<(ModificationSpecificity, f32)>,
    max_variable_mods: usize,
}

impl NonspecificDigest {
    fn new(parameters: &Parameters) -> Self {
        let enzyme: EnzymeParameters = parameters.enzyme.clone().into();
        Self {
            min_len: enzyme.min_len,
            max_len: enzyme.max_len,
            min_mass: parameters.peptide_min_mass,
            max_mass: parameters.peptide_max_mass,
            variable_mods: parameters.variable_mod_list(),
            max_variable_mods: parameters.max_variable_mods,
        }
    }

    /// `EnzymeParameters::digest` also requires a non-empty window.
    fn lengths(&self) -> RangeInclusive<usize> {
        self.min_len.max(1)..=self.max_len
    }
}

fn window_position(first: usize, length: usize, protein_length: usize) -> Position {
    match (first == 0, first + length == protein_length) {
        (true, true) => Position::Full,
        (true, false) => Position::Nterm,
        (false, true) => Position::Cterm,
        (false, false) => Position::Internal,
    }
}

/// `EnzymeParameters::digest` keeps a sequence once per protein, at its first start,
/// with that start's `Position`. A later occurrence must not yield rows when its
/// position differs from the first one's: when the first one is the protein's
/// N-terminal window, or when the later one is its C-terminal window. A later internal
/// occurrence of an internal first one yields the same rows, which are merged.
fn repeats_an_occurrence_with_another_position(
    residues: &[u8],
    first: usize,
    length: usize,
    position: Position,
) -> bool {
    if first == 0 {
        return false;
    }
    let window = &residues[first..first + length];
    if window == &residues[..length] {
        return true;
    }
    position == Position::Cterm
        && residues[1..first + length - 1]
            .windows(length)
            .any(|earlier| earlier == window)
}

/// About `count` runs of consecutive proteins holding similar residue counts.
fn protein_chunks(proteome: &Proteome, count: usize) -> Vec<Range<usize>> {
    let count = count.max(1);
    let per_chunk = ((proteome.residue_count() + count - 1) / count).max(1);
    let mut chunks = Vec::new();
    let mut chunk_start = 0;
    let mut residues = 0;
    for protein in 0..proteome.protein_count() {
        residues += proteome.protein_range(protein).len();
        if residues >= per_chunk {
            chunks.push(chunk_start..protein + 1);
            chunk_start = protein + 1;
            residues = 0;
        }
    }
    if chunk_start < proteome.protein_count() {
        chunks.push(chunk_start..proteome.protein_count());
    }
    chunks
}

/// Per-chunk scratch space for [`CompactPeptides::generate_window_rows`].
struct Scratch {
    variable_sites: Vec<(Site, f32)>,
    static_only: ModifiedForm,
    form: ModifiedForm,
}

/// A peptide's `modifications`, `nterm` and `cterm` as `Peptide` holds them, without
/// the heap. The methods replay `Peptide::apply_site`, `Peptide::static_mods` and
/// `Peptide::modification_mass`.
struct ModifiedForm {
    length: usize,
    modifications: [f32; MAX_PEPTIDE_LENGTH],
    nterm: Option<f32>,
    cterm: Option<f32>,
}

impl ModifiedForm {
    fn new() -> Self {
        Self {
            length: 0,
            modifications: [0.0; MAX_PEPTIDE_LENGTH],
            nterm: None,
            cterm: None,
        }
    }

    /// Unmodified, `length` residues long.
    fn reset(&mut self, length: usize) {
        self.length = length;
        self.modifications[..length].fill(0.0);
        self.nterm = None;
        self.cterm = None;
    }

    fn modifications(&self) -> &[f32] {
        &self.modifications[..self.length]
    }

    fn apply_site(&mut self, site: Site, mass: f32) {
        match site {
            Site::Nterm => {
                if self.nterm.is_none() {
                    self.nterm = Some(0.0 + mass);
                }
            }
            Site::Cterm => {
                if self.cterm.is_none() {
                    self.cterm = Some(0.0 + mass);
                }
            }
            Site::Sequence(index) => {
                if self.modifications[index as usize] == 0.0 {
                    self.modifications[index as usize] += mass;
                }
            }
        }
    }

    /// Every static mod, in the map's own order, as `Peptide::apply` applies them.
    fn apply_static_mods(
        &mut self,
        static_mods: &HashMap<ModificationSpecificity, f32>,
        residues: &[u8],
        position: Position,
    ) {
        for (&target, &mass) in static_mods {
            for_each_site(target, residues, position, |site| {
                self.apply_site(site, mass)
            });
        }
    }

    fn modification_mass(&self) -> f32 {
        self.modifications().iter().sum::<f32>()
            + self.nterm.unwrap_or(0.0)
            + self.cterm.unwrap_or(0.0)
    }

    /// `Peptide::reverse`'s modifications: the interior reversed, both ends kept.
    fn reverse_interior(&mut self) {
        let last = self.length.saturating_sub(1);
        if last > 1 {
            self.modifications[1..last].reverse();
        }
    }
}

/// The sites `target` takes on a window: those `Peptide::push_resi` offers a
/// variable mod and `Peptide::static_mods` fills with a static one.
fn for_each_site(
    target: ModificationSpecificity,
    residues: &[u8],
    position: Position,
    mut visit: impl FnMut(Site),
) {
    let first = *residues.first().unwrap_or(&0);
    let last = *residues.last().unwrap_or(&0);
    let last_index = residues.len().saturating_sub(1) as u32;
    let protein_nterm = matches!(position, Position::Nterm | Position::Full);
    let protein_cterm = matches!(position, Position::Cterm | Position::Full);
    match target {
        ModificationSpecificity::PeptideN(None) => visit(Site::Nterm),
        ModificationSpecificity::PeptideN(Some(residue)) if residue == first => {
            visit(Site::Sequence(0))
        }
        ModificationSpecificity::PeptideC(None) => visit(Site::Cterm),
        ModificationSpecificity::PeptideC(Some(residue)) if residue == last => {
            visit(Site::Sequence(last_index))
        }
        ModificationSpecificity::ProteinN(None) if protein_nterm => visit(Site::Nterm),
        ModificationSpecificity::ProteinN(Some(residue)) if protein_nterm && residue == first => {
            visit(Site::Sequence(0))
        }
        ModificationSpecificity::ProteinC(None) if protein_cterm => visit(Site::Cterm),
        ModificationSpecificity::ProteinC(Some(residue)) if protein_cterm && residue == last => {
            visit(Site::Sequence(last_index))
        }
        ModificationSpecificity::Residue(residue) => {
            for (index, &at) in residues.iter().enumerate() {
                if at == residue {
                    visit(Site::Sequence(index as u32));
                }
            }
        }
        _ => {}
    }
}

/// Which decoys [`Parameters::digest`] drops: a decoy whose sequence (its target's,
/// interior reversed) equals some target window. One bit per (length, start).
struct DecoyFilter {
    min_len: usize,
    residue_count: usize,
    dropped: Vec<u64>,
}

impl DecoyFilter {
    fn build(proteome: &Proteome, lengths: RangeInclusive<usize>) -> Self {
        let min_len = *lengths.start();
        let residue_count = proteome.residue_count();
        let bits = lengths.clone().count() * residue_count;
        let dropped: Vec<AtomicU64> = (0..(bits + 63) / 64)
            .into_par_iter()
            .map(|_| AtomicU64::new(0))
            .collect();
        for length in lengths {
            let targets = WindowSet::build(proteome, length);
            (0..proteome.protein_count())
                .into_par_iter()
                .for_each(|protein| {
                    let mut decoy = [0u8; MAX_PEPTIDE_LENGTH];
                    for start in window_starts(proteome, protein, length) {
                        let target = &proteome.residues[start..start + length];
                        if targets.contains(proteome, reversed_interior(target, &mut decoy)) {
                            let bit = (length - min_len) * residue_count + start;
                            dropped[bit / 64].fetch_or(1 << (bit % 64), AtomicOrdering::Relaxed);
                        }
                    }
                });
        }
        Self {
            min_len,
            residue_count,
            dropped: dropped.into_iter().map(AtomicU64::into_inner).collect(),
        }
    }

    fn is_dropped(&self, start: usize, length: usize) -> bool {
        let bit = (length - self.min_len) * self.residue_count + start;
        (self.dropped[bit / 64] >> (bit % 64)) & 1 == 1
    }
}

/// `Peptide::reverse`'s sequence of `target`, written into `buffer`.
fn reversed_interior<'a>(target: &[u8], buffer: &'a mut [u8; MAX_PEPTIDE_LENGTH]) -> &'a [u8] {
    let decoy = &mut buffer[..target.len()];
    decoy.copy_from_slice(target);
    let last = target.len().saturating_sub(1);
    if last > 1 {
        decoy[1..last].reverse();
    }
    decoy
}

/// Proteome offsets of `protein`'s windows of `length` residues.
fn window_starts(proteome: &Proteome, protein: usize, length: usize) -> Range<usize> {
    let range = proteome.protein_range(protein);
    match range.len() >= length {
        true => range.start..range.end - length + 1,
        false => range.start..range.start,
    }
}

/// Target windows of one length: an open-addressing set of their proteome offsets,
/// keyed by their residues.
struct WindowSet {
    length: usize,
    shift: u32,
    slots: Vec<AtomicU32>,
}

impl WindowSet {
    const EMPTY: u32 = u32::MAX;

    fn build(proteome: &Proteome, length: usize) -> Self {
        let windows: usize = (0..proteome.protein_count())
            .map(|protein| window_starts(proteome, protein, length).len())
            .sum();
        let capacity = (2 * windows).next_power_of_two().max(2);
        let set = Self {
            length,
            shift: 64 - capacity.trailing_zeros(),
            slots: (0..capacity)
                .into_par_iter()
                .map(|_| AtomicU32::new(Self::EMPTY))
                .collect(),
        };
        (0..proteome.protein_count())
            .into_par_iter()
            .for_each(|protein| {
                for start in window_starts(proteome, protein, length) {
                    set.insert(proteome, start as u32);
                }
            });
        set
    }

    fn first_slot(&self, residues: &[u8]) -> usize {
        let mut hasher = FnvHasher::default();
        hasher.write(residues);
        (hasher.finish().wrapping_mul(0x9E37_79B9_7F4A_7C15) >> self.shift) as usize
    }

    fn insert(&self, proteome: &Proteome, start: u32) {
        let window = proteome.window(start, self.length as u8);
        let mask = self.slots.len() - 1;
        let mut slot = self.first_slot(window);
        loop {
            let occupant = match self.slots[slot].compare_exchange(
                Self::EMPTY,
                start,
                AtomicOrdering::AcqRel,
                AtomicOrdering::Acquire,
            ) {
                Ok(_) => return,
                Err(occupant) => occupant,
            };
            if proteome.window(occupant, self.length as u8) == window {
                return;
            }
            slot = (slot + 1) & mask;
        }
    }

    fn contains(&self, proteome: &Proteome, residues: &[u8]) -> bool {
        let mask = self.slots.len() - 1;
        let mut slot = self.first_slot(residues);
        loop {
            let occupant = self.slots[slot].load(AtomicOrdering::Relaxed);
            if occupant == Self::EMPTY {
                return false;
            }
            if proteome.window(occupant, self.length as u8) == residues {
                return true;
            }
            slot = (slot + 1) & mask;
        }
    }
}

/// Fine mass bins over the digest's mass range. Equal masses share a bin, so a
/// group of duplicates never straddles two.
struct MassBins {
    min: f64,
    width: f64,
    count: usize,
}

impl MassBins {
    fn new(min_mass: f32, max_mass: f32) -> Self {
        let span = (max_mass as f64 - min_mass as f64).max(0.0);
        let width = FINE_BIN_DA.max(span / MAX_FINE_BINS as f64);
        Self {
            min: min_mass as f64,
            width,
            count: (span / width) as usize + 1,
        }
    }

    fn of(&self, mass: f32) -> usize {
        (((mass as f64 - self.min) / self.width) as usize).min(self.count - 1)
    }
}

/// Consecutive runs of fine bins holding about [`ROWS_PER_PLACEMENT_BIN`] rows each,
/// laid out in mass order in the row buffer.
struct Placement {
    bin_of_fine: Vec<u32>,
    starts: Vec<usize>,
}

impl Placement {
    fn new(fine_counts: &[u64]) -> Self {
        let mut bin_of_fine = Vec::with_capacity(fine_counts.len());
        let mut starts = vec![0];
        let mut total = 0;
        let mut in_bin = 0;
        for &count in fine_counts {
            if in_bin >= ROWS_PER_PLACEMENT_BIN && count > 0 {
                starts.push(total);
                in_bin = 0;
            }
            bin_of_fine.push((starts.len() - 1) as u32);
            in_bin += count;
            total += count as usize;
        }
        starts.push(total);
        Self {
            bin_of_fine,
            starts,
        }
    }

    fn bins(&self) -> usize {
        self.starts.len() - 1
    }

    fn total(&self) -> usize {
        self.starts[self.bins()]
    }

    fn range(&self, bin: usize) -> Range<usize> {
        self.starts[bin]..self.starts[bin + 1]
    }
}

/// A row buffer many threads write into, each at indices it alone claimed.
struct RowSlots<'a> {
    rows: *mut CompactPeptide,
    len: usize,
    _rows: PhantomData<&'a mut [CompactPeptide]>,
}

// SAFETY: threads only write, each to indices no other thread touches (see `write`).
unsafe impl Sync for RowSlots<'_> {}

impl<'a> RowSlots<'a> {
    fn new(rows: &'a mut [CompactPeptide]) -> Self {
        Self {
            rows: rows.as_mut_ptr(),
            len: rows.len(),
            _rows: PhantomData,
        }
    }

    /// # Safety
    /// No other thread reads or writes `index` while the buffer is shared.
    unsafe fn write(&self, index: usize, row: CompactPeptide) {
        assert!(index < self.len);
        self.rows.add(index).write(row);
    }
}
