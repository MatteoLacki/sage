use crate::enzyme::{group_digests, Enzyme, EnzymeParameters};
use crate::fasta::Fasta;
use crate::ion_series::{IonSeries, Kind};
use crate::mass::Tolerance;
use crate::modification::{validate_mods, validate_var_mods, ModificationSpecificity};
use crate::peptide::Peptide;
use dashmap::DashSet;
use fnv::FnvBuildHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::cell::RefCell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::hash::Hash;

#[derive(Deserialize, Serialize, Clone, Debug)]
pub struct EnzymeBuilder {
    /// How many missed cleavages to use
    pub missed_cleavages: Option<u8>,
    /// Minimum peptide length that will be fragmented
    pub min_len: Option<usize>,
    /// Maximum peptide length that will be fragmented
    pub max_len: Option<usize>,
    pub cleave_at: Option<String>,
    pub restrict: Option<String>,
    pub c_terminal: Option<bool>,
    pub semi_enzymatic: Option<bool>,
}

impl Default for EnzymeBuilder {
    fn default() -> Self {
        Self {
            missed_cleavages: Some(0),
            min_len: Some(5),
            max_len: Some(50),
            cleave_at: Some("KR".into()),
            restrict: Some("P".into()),
            c_terminal: Some(true),
            semi_enzymatic: Some(false),
        }
    }
}

impl From<EnzymeBuilder> for EnzymeParameters {
    fn from(en: EnzymeBuilder) -> EnzymeParameters {
        EnzymeParameters {
            missed_cleavages: en.missed_cleavages.unwrap_or(1),
            min_len: en.min_len.unwrap_or(5),
            max_len: en.max_len.unwrap_or(50),
            enzyme: Enzyme::new(
                &en.cleave_at.unwrap_or_else(|| "KR".into()),
                &en.restrict.unwrap_or_else(|| "".into()),
                en.c_terminal.unwrap_or(true),
                en.semi_enzymatic.unwrap_or(false),
            ),
        }
    }
}

#[derive(Deserialize, Default)]
/// Parameters used for generating the fragment database
pub struct Builder {
    /// This parameter allows tuning of the internal search structure
    pub bucket_size: Option<usize>,

    pub enzyme: Option<EnzymeBuilder>,
    /// Minimum peptide monoisotopic mass that will be fragmented
    pub peptide_min_mass: Option<f32>,
    /// Maximum peptide monoisotopic mass that will be fragmented
    pub peptide_max_mass: Option<f32>,
    /// Which kind of fragment ions to generate (a, b, c, x, y, z)
    pub ion_kinds: Option<Vec<Kind>>,
    /// Minimum ion index to be generated: 1 will remove b1/y1 ions
    /// 2 will remove b1/b2/y1/y2 ions, etc
    pub min_ion_index: Option<usize>,
    /// Static modifications to add to matching amino acids
    pub static_mods: Option<HashMap<String, f32>>,
    /// Variable modifications to add to matching amino acids
    pub variable_mods: Option<HashMap<String, Vec<f32>>>,
    /// Limit number of variable modifications on a peptide
    pub max_variable_mods: Option<usize>,
    /// Use this prefix for decoy proteins
    pub decoy_tag: Option<String>,

    pub generate_decoys: Option<bool>,
    /// Path to fasta database
    pub fasta: Option<String>,
    /// Number of sequences to handle simultaneously when pre-filtering the db
    pub prefilter_chunk_size: Option<usize>,
    /// Pre-filter the database to minimize memory usage
    pub prefilter: Option<bool>,
    /// Pre-filter the database with a minimal amount of memory at the cost of speed
    pub prefilter_low_memory: Option<bool>,
}

impl Builder {
    pub fn make_parameters(self) -> Parameters {
        let bucket_size = self.bucket_size.unwrap_or(8192).next_power_of_two();
        Parameters {
            bucket_size,
            peptide_min_mass: self.peptide_min_mass.unwrap_or(500.0),
            peptide_max_mass: self.peptide_max_mass.unwrap_or(5000.0),
            ion_kinds: self.ion_kinds.unwrap_or(vec![Kind::B, Kind::Y]),
            min_ion_index: self.min_ion_index.unwrap_or(2),
            decoy_tag: self.decoy_tag.unwrap_or_else(|| "rev_".into()),
            enzyme: self.enzyme.unwrap_or_default(),
            static_mods: validate_mods(self.static_mods),
            variable_mods: validate_var_mods(self.variable_mods),
            max_variable_mods: self.max_variable_mods.map(|x| x.max(1)).unwrap_or(2),
            generate_decoys: self.generate_decoys.unwrap_or(true),
            fasta: self.fasta.expect("A fasta file must be provided!"),
            prefilter_chunk_size: self.prefilter_chunk_size.unwrap_or(0),
            prefilter: self.prefilter.unwrap_or(false),
            prefilter_low_memory: self.prefilter_low_memory.unwrap_or(true),
        }
    }

    pub fn update_fasta(&mut self, fasta: String) {
        self.fasta = Some(fasta)
    }
}

#[derive(Serialize, Clone, Debug)]
pub struct Parameters {
    pub bucket_size: usize,
    pub enzyme: EnzymeBuilder,
    pub peptide_min_mass: f32,
    pub peptide_max_mass: f32,
    pub ion_kinds: Vec<Kind>,
    pub min_ion_index: usize,
    pub static_mods: HashMap<ModificationSpecificity, f32>,
    pub variable_mods: HashMap<ModificationSpecificity, Vec<f32>>,
    pub max_variable_mods: usize,
    pub decoy_tag: String,
    pub generate_decoys: bool,
    pub fasta: String,
    pub prefilter_chunk_size: usize,
    pub prefilter: bool,
    pub prefilter_low_memory: bool,
}

impl Parameters {
    pub fn auto_calculate_prefilter_chunk_size(&mut self, fasta: &Fasta) {
        const MAX_PEPS_PER_CHUNK: usize = 2usize.pow(23);
        self.prefilter_chunk_size = match self.prefilter_chunk_size {
            0 => {
                let enzyme = self.enzyme.clone().into();
                let total_unmodified_pep_count: usize = fasta.digest(&enzyme).len();
                let mod_count_estimate =
                    (self.variable_mods.len() + 1) * (1 << self.max_variable_mods);
                let chunk_count =
                    mod_count_estimate * total_unmodified_pep_count / MAX_PEPS_PER_CHUNK;
                if chunk_count == 0 {
                    fasta.targets.len()
                } else {
                    fasta.targets.len() / chunk_count
                }
            }
            x => x,
        };
    }

    pub fn digest(&self, fasta: &Fasta) -> Vec<Peptide> {
        log::trace!("digesting fasta");
        let enzyme = self.enzyme.clone().into();
        // Generate all tryptic peptide sequences, including reversed (decoy)
        // and missed cleavages, if applicable.
        let digests = fasta.digest(&enzyme);

        log::trace!("grouping digests");
        let start_num = digests.len();
        let digests = group_digests(digests);
        log::trace!(
            "grouped {} digests into {} groups",
            start_num,
            digests.len()
        );

        let mods = self
            .variable_mods
            .iter()
            .flat_map(|(a, b)| b.iter().map(|b| (*a, *b)))
            .collect::<Vec<_>>();

        let targets: DashSet<_, FnvBuildHasher> = DashSet::default();
        digests
            .par_iter()
            .filter(|digest| !digest.reference.decoy)
            .for_each(|digest| {
                targets.insert(digest.reference.sequence.clone().into_bytes());
            });

        log::trace!("modifying peptides");
        let mut target_decoys = digests
            .into_par_iter()
            .map(Peptide::try_from)
            .filter_map(Result::ok)
            .flat_map_iter(|peptide| {
                peptide
                    .apply(&mods, &self.static_mods, self.max_variable_mods)
                    .into_iter()
                    .filter(|peptide| {
                        peptide.monoisotopic >= self.peptide_min_mass
                            && peptide.monoisotopic <= self.peptide_max_mass
                    })
                    .flat_map(|peptide| {
                        if self.generate_decoys {
                            vec![peptide.reverse(), peptide].into_iter()
                        } else {
                            vec![peptide].into_iter()
                        }
                    })
                    .filter(|peptide| !peptide.decoy || !targets.contains(&(peptide.sequence[..])))
            })
            .collect::<Vec<_>>();

        Self::reorder_peptides(&mut target_decoys);

        target_decoys
    }

    pub fn reorder_peptides(target_decoys: &mut Vec<Peptide>) {
        log::trace!("sorting and deduplicating peptides");

        let init_size = target_decoys.len();
        // This is equivalent to a stable sort
        target_decoys.par_sort_unstable_by(|a, b| {
            a.monoisotopic
                .total_cmp(&b.monoisotopic)
                .then_with(|| a.initial_sort(b))
        });
        target_decoys.dedup_by(|remove, keep| {
            if remove.monoisotopic == keep.monoisotopic
                && remove.sequence == keep.sequence
                && remove.modifications == keep.modifications
                && remove.nterm == keep.nterm
                && remove.cterm == keep.cterm
            {
                keep.proteins.extend(remove.proteins.iter().cloned());
                // When merging peptides from different Fastas,
                // decoys in one fasta might be targets in another
                keep.decoy &= remove.decoy;
                true
            } else {
                false
            }
        });

        target_decoys
            .par_iter_mut()
            .for_each(|peptide| peptide.proteins.sort_unstable());

        let num_dropped = init_size - target_decoys.len();
        log::trace!(
            "dropped {} t/d pairs, remaining {}",
            num_dropped,
            target_decoys.len(),
        );
    }

    pub fn build(self, fasta: Fasta) -> IndexedDatabase {
        let target_decoys = self.digest(&fasta);
        self.build_from_peptides(target_decoys)
    }

    pub fn build_from_peptides(self, target_decoys: Vec<Peptide>) -> IndexedDatabase {
        log::trace!("generating fragments");

        // Finally, perform in silico digest for our target sequences
        // Note that multiple charge states are actually handled by
        // [`SpectrumProcessor`] or during scoring - all theoretical
        // fragments are monoisotopic/uncharged
        let mut fragments = target_decoys
            .par_iter()
            .enumerate()
            .flat_map_iter(|(idx, peptide)| {
                // Generate both B and Y ions, then filter down to make sure that
                // theoretical fragments are within the search space
                self.ion_kinds
                    .iter()
                    .flat_map(|kind| IonSeries::new(peptide, *kind).enumerate())
                    .filter(|(ion_idx, ion)| {
                        // Don't store b1, b2, y1, y2 ions for preliminary scoring

                        match ion.kind {
                            Kind::A | Kind::B | Kind::C => (ion_idx + 1) > self.min_ion_index,
                            Kind::X | Kind::Y | Kind::Z => {
                                peptide.sequence.len().saturating_sub(1) - ion_idx
                                    > self.min_ion_index
                            }
                        }
                    })
                    .map(move |(_, ion)| Theoretical {
                        peptide_index: PeptideIx(idx as u32),
                        fragment_mz: ion.monoisotopic_mass,
                    })
            })
            .collect::<Vec<_>>();
        log::trace!("finalizing index");

        // Sort all of our theoretical fragments by m/z, from low to high
        fragments.par_sort_unstable_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));

        // Now, we bucket all of our theoretical fragments, and within each bucket
        // sort by precursor m/z - and save the minimum *fragment* m/z in a separate
        // vector so that we can perform an efficient binary search to reduce
        // the number of in silico fragments we evaluate
        //
        // Imagine our theoretical fragments look like this
        //
        // Fragment        A      B       C       D       E       F       G       H
        // Fragment m/z [ 1.0    1.2     1.3     2.5     2.5     2.6     3.5     4.0 ]
        // Parent m/z   [ 500    439     291     800     142     515     517     232 ]
        //
        // If we apply a bucket size of 4 we will end up with the following:
        //
        // Fragment        C      B       A       D       E       H       F       G
        // Fragment m/z [ 1.3    1.2     1.0     2.5     2.5     4.0     3.5     2.6 ]
        // Parent m/z   [ 291    439     500     800     142     232     515     517 ]
        //              |___________________________|   |____________________________|
        //               Bucket 1: min m/z 1.0          Bucket 2: min m/z 2.5
        //
        // * Example query: Fragment m/z 1.3 - 1.9 & Precursor m/z: 450 - 900
        // 1) Perform a binary search to narrow down our window to Bucket 1 only
        //      * Bucket 2 has a min m/z outside of our query range - nothing here can match
        //
        // Fragment        C      B       A       D
        // Fragment m/z [ 1.3    1.2     1.0     2.5
        // Parent m/z   [ 291    439     500     800
        //                            |_____________|
        //                                    ^
        //                                    |
        // Window with matching precursors ___|

        // and within Bucket 1, we can perform another binary search to find fragments
        // matching our desired precursor m/z tolerance

        let min_value = fragments
            .par_chunks_mut(self.bucket_size)
            .map(|chunk| {
                // There should always be at least one item in the chunk!
                //  we know the chunk is already sorted by fragment_mz too, so this is minimum value
                let min = chunk[0].fragment_mz;
                chunk.par_sort_unstable_by(|a, b| a.peptide_index.cmp(&b.peptide_index));
                min
            })
            .collect::<Vec<_>>();

        let potential_mods = self
            .variable_mods
            .iter()
            .flat_map(|(a, b)| b.iter().map(|b| (*a, *b)))
            .collect::<Vec<(ModificationSpecificity, f32)>>();

        let (fragment_mzs, fragment_peptide_indices) = fragments
            .into_par_iter()
            .map(|fragment| (fragment.fragment_mz, fragment.peptide_index))
            .unzip();

        IndexedDatabase {
            peptides: target_decoys,
            fragment_mzs,
            fragment_peptide_indices,
            min_value,
            bucket_size: self.bucket_size,
            ion_kinds: self.ion_kinds,
            generate_decoys: self.generate_decoys,
            potential_mods,
            decoy_tag: self.decoy_tag,
        }
    }
}

#[derive(Hash, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize)]
#[repr(transparent)]
pub struct PeptideIx(pub u32);

// This is unsafe for use outside of this crate
impl Default for PeptideIx {
    fn default() -> Self {
        Self(u32::MAX)
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize)]
pub struct Theoretical {
    pub peptide_index: PeptideIx,
    pub fragment_mz: f32,
}

#[derive(Default)]
pub struct IndexedDatabase {
    pub peptides: Vec<Peptide>,
    /// Parallel columns, ordered by peptide index within each fragment-mass page.
    pub fragment_mzs: Vec<f32>,
    pub fragment_peptide_indices: Vec<PeptideIx>,
    pub ion_kinds: Vec<Kind>,
    pub min_value: Vec<f32>,
    /// Keep a list of potential (AA, mass) modifications for RT prediction
    pub potential_mods: Vec<(ModificationSpecificity, f32)>,
    pub bucket_size: usize,
    pub generate_decoys: bool,
    pub decoy_tag: String,
}

impl IndexedDatabase {
    /// Create a new [`IndexedQuery`] for a specific [`ProcessedSpectrum`]
    ///
    /// All matches returned by the query will be within the specified tolerance
    /// parameters
    pub fn query(&self, precursor_mass: f32, precursor_tol: Tolerance) -> IndexedQuery<'_> {
        let (precursor_lo, precursor_hi) = precursor_tol.bounds(precursor_mass);

        let (pre_idx_lo, pre_idx_hi) = binary_search_slice(
            &self.peptides,
            |p, bounds| p.monoisotopic.total_cmp(bounds),
            precursor_lo,
            precursor_hi,
        );

        IndexedQuery {
            db: self,
            precursor_mass,
            precursor_tol,
            pre_idx_lo,
            pre_idx_hi,
        }
    }

    pub fn size(&self) -> usize {
        self.fragment_mzs.len()
    }

    pub fn fragments(&self) -> impl Iterator<Item = Theoretical> + '_ {
        self.fragment_mzs
            .iter()
            .zip(&self.fragment_peptide_indices)
            .map(|(&fragment_mz, &peptide_index)| Theoretical {
                peptide_index,
                fragment_mz,
            })
    }

    pub fn buckets(&self) -> &[f32] {
        &self.min_value
    }

    pub fn serialize(&self) {
        use std::io::Write;
        let mut wtr = std::io::BufWriter::new(std::fs::File::create("fragments.bin").unwrap());
        for fragment in self.fragments() {
            let _ = wtr.write(&fragment.fragment_mz.to_le_bytes()).unwrap();
            let _ = wtr.write(&fragment.peptide_index.0.to_le_bytes()).unwrap();
        }
        wtr.flush().unwrap();

        let mut wtr = std::io::BufWriter::new(std::fs::File::create("peptides.csv").unwrap());
        writeln!(wtr, "peptide,proteins,monoisotopic,decoy").unwrap();
        for fragment in &self.peptides {
            writeln!(
                wtr,
                "{},{},{},{}",
                fragment,
                fragment.proteins(&self.decoy_tag, self.generate_decoys),
                fragment.monoisotopic,
                fragment.decoy
            )
            .unwrap();
        }
        wtr.flush().unwrap();
    }
}

impl std::ops::Index<PeptideIx> for IndexedDatabase {
    type Output = Peptide;

    fn index(&self, index: PeptideIx) -> &Self::Output {
        &self.peptides[index.0 as usize]
    }
}

pub struct IndexedQuery<'d> {
    db: &'d IndexedDatabase,
    precursor_mass: f32,
    precursor_tol: Tolerance,
    pub pre_idx_lo: usize,
    pub pre_idx_hi: usize,
}

impl IndexedQuery<'_> {
    /// Search for a specified `fragment_mz` within the database, using an
    /// explicit `fragment_tol` (rather than a fixed one baked into the
    /// query) so callers can vary it per call.
    pub fn page_search(
        &self,
        mass: f32,
        fragment_tol: Tolerance,
    ) -> impl Iterator<Item = Theoretical> + '_ {
        let (fragment_lo, fragment_hi) = fragment_tol.bounds(mass);
        let (precursor_lo, precursor_hi) = self.precursor_tol.bounds(self.precursor_mass);

        // Locate the left and right page indices that contain matching fragments
        // Note that we need to multiply by `bucket_size` to transform these into
        // indices into the parallel fragment columns
        let (left_idx, right_idx) = binary_search_slice(
            &self.db.min_value,
            |min, bounds| min.total_cmp(bounds),
            fragment_lo,
            fragment_hi,
        );

        // It is absolutely critical that we do not cross page boundaries!
        // If we do, we can no longer rely on total ordering of peptide_index (precursor m/z)
        (left_idx..right_idx).flat_map(move |page| {
            let left_idx = page * self.db.bucket_size;
            // Last chunk not guaranted to be modulo bucket size, make sure we don't
            // accidentally go out of bounds!
            let right_idx = ((page + 1) * self.db.bucket_size).min(self.db.size());

            // Narrow down into our region of interest, then perform another binary
            // search to further refine down to the slice of matching precursor mzs
            let slice = &self.db.fragment_peptide_indices[left_idx..right_idx];

            let (inner_left, inner_right) = binary_search_slice(
                slice,
                |idx, bounds| (idx.0 as usize).cmp(bounds),
                self.pre_idx_lo,
                self.pre_idx_hi,
            );

            // Finally, filter down our slice into exact matches only
            (inner_left..inner_right)
                .map(move |i| Theoretical {
                    peptide_index: slice[i],
                    fragment_mz: self.db.fragment_mzs[left_idx + i],
                })
                .filter(move |frag| {
                    // This looks somewhat complicated, but it's a consequence of
                    // how the `binary_search_slice` function works - it will return
                    // the set of indices that maximally cover the desired range - the exact
                    // `left` and `right` indices may be valid, or just outside of the range.
                    // Anything interior of `left` and `right` is guaranteed to be within the
                    // precursor tolerance, so we just need to check the edge cases
                    //
                    // Previously, a direct lookup to check the mass of the current fragment was
                    // performed, but the pointer indirection + float comparison can slow down
                    // open searches by as much as 2x!!
                    // e.g. used to be `self.db[frag.peptide_index].monoisotopic >= precursor_lo`
                    (frag.peptide_index.0 > self.pre_idx_lo as u32
                        || (frag.peptide_index.0 == self.pre_idx_lo as u32
                            && self.db[frag.peptide_index].monoisotopic >= precursor_lo))
                        && (frag.peptide_index.0 < self.pre_idx_hi as u32
                            || (frag.peptide_index.0 == self.pre_idx_hi as u32
                                && self.db[frag.peptide_index].monoisotopic <= precursor_hi))
                        && frag.fragment_mz >= fragment_lo
                        && frag.fragment_mz <= fragment_hi
                })
        })
    }

    /// Batched `page_search`: search many `(mass, fragment_tol)` windows
    /// (e.g. one spectrum's whole peak list x fragment-charge range) in one
    /// call, computing each page's precursor-scoped `[inner_left,
    /// inner_right)` binary search only once and reusing it across every
    /// window that overlaps that page -- repeated `page_search` calls
    /// recompute this independently per window even though it only depends
    /// on `(page, pre_idx_lo, pre_idx_hi)`, and `pre_idx_lo`/`pre_idx_hi`
    /// are fixed for this whole `IndexedQuery`. Windows need not be
    /// pre-sorted (page numbers are sorted internally), though in practice
    /// callers already have them sorted by mass (spectra peaks are), which
    /// keeps that sort near-free.
    ///
    /// Takes a callback instead of returning a `Vec` -- matches stream to
    /// `on_match` as found, same zero-allocation-for-results shape
    /// `page_search`'s lazy iterator chain had (a first version collected
    /// into an owned `Vec<&Theoretical>`, which turned out to cost more,
    /// on typical low-peak-sharing spectra, than the inner-search reuse
    /// saved -- see `docs/ai/reuse_index_bins.md`). The two remaining
    /// internal scratch buffers (`bounds`, `page_window`) are thread-local
    /// and reused across calls -- rayon's worker threads are long-lived
    /// OS threads, so this persists correctly across the many spectra one
    /// worker processes, without needing to thread a buffer through every
    /// call site by hand.
    pub fn page_search_batch(
        &self,
        windows: &[(f32, Tolerance)],
        mut on_match: impl FnMut(Theoretical),
    ) {
        let (precursor_lo, _) = self.precursor_tol.bounds(self.precursor_mass);
        // binary_search_slice includes one predecessor at the lower boundary.
        // Convert to an exact half-open peptide-ID range before scanning pages.
        let exact_pre_idx_lo = self.pre_idx_lo
            + usize::from(
                self.db
                    .peptides
                    .get(self.pre_idx_lo)
                    .map_or(false, |p| p.monoisotopic < precursor_lo),
            );
        #[cfg(target_arch = "x86_64")]
        let use_avx2 = std::arch::is_x86_feature_detected!("avx2");

        PAGE_SEARCH_SCRATCH.with(|scratch| {
            let mut scratch = scratch.borrow_mut();
            let PageSearchScratch {
                bounds,
                page_window,
                groups,
                page_bounds,
            } = &mut *scratch;
            bounds.clear();
            page_window.clear();

            for (wi, &(mass, fragment_tol)) in windows.iter().enumerate() {
                let (fragment_lo, fragment_hi) = fragment_tol.bounds(mass);
                bounds.push((fragment_lo, fragment_hi));

                let (left_page, right_page) = binary_search_slice(
                    &self.db.min_value,
                    |min, b| min.total_cmp(b),
                    fragment_lo,
                    fragment_hi,
                );
                for page in left_page..right_page {
                    page_window.push((page, wi));
                }
            }
            page_window.sort_unstable_by_key(|&(page, _)| page);

            groups.clear();
            let mut i = 0;
            while i < page_window.len() {
                let page = page_window[i].0;
                let mut j = i;
                while j < page_window.len() && page_window[j].0 == page {
                    j += 1;
                }
                groups.push((page, i, j));
                i = j;
            }

            resolve_page_peptide_bounds(
                &self.db.fragment_peptide_indices,
                self.db.bucket_size,
                self.db.size(),
                groups,
                exact_pre_idx_lo,
                self.pre_idx_hi,
                page_bounds,
            );

            for (group_idx, &(page, window_start, window_end)) in groups.iter().enumerate() {
                let (inner_left, inner_right) = page_bounds[group_idx];
                let (inner_left, inner_right) = (inner_left as usize, inner_right as usize);
                let left_idx = page * self.db.bucket_size;
                let masses = &self.db.fragment_mzs[left_idx + inner_left..left_idx + inner_right];
                let ids = &self.db.fragment_peptide_indices
                    [left_idx + inner_left..left_idx + inner_right];

                let mut j = window_start;
                while j < window_end {
                    let (fragment_lo, fragment_hi) = bounds[page_window[j].1];
                    let mut emit = |idx: usize| {
                        on_match(Theoretical {
                            peptide_index: ids[idx],
                            fragment_mz: masses[idx],
                        })
                    };
                    #[cfg(target_arch = "x86_64")]
                    if use_avx2 && masses.len() >= 8 {
                        // SAFETY: runtime feature detection above guards this call.
                        unsafe { scan_masses_avx2(masses, fragment_lo, fragment_hi, &mut emit) };
                        j += 1;
                        continue;
                    }
                    scan_masses_scalar(masses, fragment_lo, fragment_hi, &mut emit);
                    j += 1;
                }
            }
        });
    }
}

/// Resolve each page's precursor-scoped `[inner_left, inner_right)` peptide-ID
/// range, running many pages' binary searches interleaved rather than one
/// after another.
///
/// Each search is ~15 *dependent* loads into a DRAM-resident ID column
/// (198M fragments here), so a search run on its own stalls on memory for
/// nearly its whole duration and leaves the core's miss parallelism idle.
/// Stepping `LANES` independent searches in lockstep keeps that many misses
/// in flight at once. Measured on F9477: this resolution was 94.1% of search
/// cycles, at ~200 cycles per binary-search step -- i.e. full memory latency,
/// unoverlapped. See `docs/ai/simd.md`.
///
/// `pre_idx_lo`/`pre_idx_hi` are fixed for one `IndexedQuery`, so every page
/// searches for the same two keys. `inner_right` is resolved against the
/// whole page instead of `ids[inner_left..]`: the column is sorted ascending
/// and `pre_idx_lo <= pre_idx_hi`, so both give the same index, and searching
/// the full page makes the two searches independent -- doubling how many can
/// be in flight.
fn resolve_page_peptide_bounds(
    all_ids: &[PeptideIx],
    bucket_size: usize,
    total: usize,
    groups: &[(usize, usize, usize)],
    pre_idx_lo: usize,
    pre_idx_hi: usize,
    out: &mut Vec<(u32, u32)>,
) {
    /// Pages resolved per chunk; each contributes two independent searches,
    /// so `LANES * 2` misses are in flight at once. Sized against the core's
    /// outstanding-miss capacity, not the vector width.
    const LANES: usize = 16;
    const SLOTS: usize = LANES * 2;

    out.clear();
    out.reserve(groups.len());
    if all_ids.is_empty() {
        return;
    }

    for chunk in groups.chunks(LANES) {
        // Two slots per page: even = lower bound, odd = upper bound.
        let mut base = [0usize; SLOTS];
        let mut remaining = [0usize; SLOTS];
        let mut page_start = [0usize; LANES];
        let mut longest_page = 0usize;

        for (c, &(page, _, _)) in chunk.iter().enumerate() {
            let left = page * bucket_size;
            let right = ((page + 1) * bucket_size).min(total);
            let len = right - left;
            page_start[c] = left;
            base[2 * c] = left;
            base[2 * c + 1] = left;
            remaining[2 * c] = len;
            remaining[2 * c + 1] = len;
            longest_page = longest_page.max(len);
        }

        // Every slot takes the same number of halving steps, so the trip
        // count is fixed and the body unrolls: a branch per slot per step
        // would serialize exactly the loads this is here to overlap. A slot
        // that has already converged (`remaining <= 1`, including the unused
        // slots of a short final chunk) gets `half == 0`, which re-probes its
        // own `base` and applies no movement.
        let steps = usize::BITS - longest_page.saturating_sub(1).leading_zeros();
        for _ in 0..steps {
            for slot in 0..SLOTS {
                let half = remaining[slot] / 2;
                let key = if slot % 2 == 0 { pre_idx_lo } else { pre_idx_hi };
                // SAFETY: `base[slot]` indexes inside its page and
                // `half <= remaining[slot]`, so `base + half - 1` stays
                // within that page; converged slots re-read `base[slot]`,
                // which is itself in bounds. `all_ids` is non-empty.
                let probe =
                    unsafe { all_ids.get_unchecked(base[slot] + half.saturating_sub(1)).0 } as usize;
                base[slot] += ((probe < key) as usize) * half;
                remaining[slot] -= half;
            }
        }

        for c in 0..chunk.len() {
            let lo_base = base[2 * c];
            let hi_base = base[2 * c + 1];
            let inner_left =
                lo_base + ((all_ids[lo_base].0 as usize) < pre_idx_lo) as usize - page_start[c];
            let inner_right =
                hi_base + ((all_ids[hi_base].0 as usize) < pre_idx_hi) as usize - page_start[c];
            out.push((inner_left as u32, inner_right as u32));
        }
    }
}

#[inline]
fn scan_masses_scalar(masses: &[f32], lo: f32, hi: f32, mut on_match: impl FnMut(usize)) {
    for (i, &mass) in masses.iter().enumerate() {
        if mass >= lo && mass <= hi {
            on_match(i);
        }
    }
}

#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn scan_masses_avx2(masses: &[f32], lo: f32, hi: f32, mut on_match: impl FnMut(usize)) {
    use std::arch::x86_64::*;

    let lower = _mm256_set1_ps(lo);
    let upper = _mm256_set1_ps(hi);
    let mut i = 0;
    while i + 8 <= masses.len() {
        // SAFETY: eight elements remain; loadu does not require alignment.
        let values = _mm256_loadu_ps(masses.as_ptr().add(i));
        let accepted = _mm256_and_ps(
            _mm256_cmp_ps(values, lower, _CMP_GE_OQ),
            _mm256_cmp_ps(values, upper, _CMP_LE_OQ),
        );
        let mut mask = _mm256_movemask_ps(accepted) as u32;
        // Emit in index order, including repeated peptide IDs in different lanes.
        while mask != 0 {
            on_match(i + mask.trailing_zeros() as usize);
            mask &= mask - 1;
        }
        i += 8;
    }
    scan_masses_scalar(&masses[i..], lo, hi, |offset| on_match(i + offset));
}

#[derive(Default)]
struct PageSearchScratch {
    bounds: Vec<(f32, f32)>,
    page_window: Vec<(usize, usize)>,
    /// One entry per distinct page: `(page, window_start, window_end)`.
    groups: Vec<(usize, usize, usize)>,
    /// `(inner_left, inner_right)` per entry in `groups`.
    page_bounds: Vec<(u32, u32)>,
}

thread_local! {
    static PAGE_SEARCH_SCRATCH: RefCell<PageSearchScratch> = RefCell::new(PageSearchScratch::default());
}

/// Return the widest `left` and `right` indices into a `slice` (sorted by the
/// function `key`) such that all values between `low` and `high` are
/// contained in `slice[left..right]`
///
/// # Invariants
///
/// * `slice[left] <= low || left == 0`
/// * `slice[right] > high || right == slice.len()`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = slice
        .partition_point(|a| key(a, &low) == Ordering::Less)
        .saturating_sub(1);

    let right_idx =
        slice[left_idx..].partition_point(|a| key(a, &high) != Ordering::Greater) + left_idx;

    (left_idx, right_idx)
}

#[cfg(test)]
mod test {
    use std::sync::Arc;

    use super::*;

    #[test]
    fn batched_pages_match_brute_force_with_boundaries_and_duplicate_hits() {
        let peptides: Vec<_> = [500.0, 500.0, 600.0, 600.0, 800.0, 900.0]
            .into_iter()
            .map(|monoisotopic| Peptide {
                monoisotopic,
                ..Default::default()
            })
            .collect();
        let mut fragments = Vec::new();
        for ordinal in 0..32 {
            for peptide in 0..peptides.len() {
                fragments.push(Theoretical {
                    peptide_index: PeptideIx(peptide as u32),
                    fragment_mz: 100.0 + ordinal as f32 * 5.0 + peptide as f32 * 0.25,
                });
            }
        }
        // Duplicate fragments must increment the same peptide more than once.
        fragments.push(fragments[0]);
        fragments.sort_by(|a, b| a.fragment_mz.total_cmp(&b.fragment_mz));
        let windows = [
            (100.0, Tolerance::Da(0.0, 0.0)),
            (120.0, Tolerance::Da(-20.0, 20.0)),
            (120.0, Tolerance::Da(-20.0, 20.0)), // repeated window retains multiplicity
            (200.0, Tolerance::Ppm(-20.0, 20.0)),
            (1000.0, Tolerance::Da(-2000.0, 2000.0)),
            (3000.0, Tolerance::Da(0.0, 1.0)),
        ];
        for bucket_size in [1, 7, 8, 9, 31, 64, 256] {
            let mut ordered = fragments.clone();
            let mut min_value = Vec::new();
            for page in ordered.chunks_mut(bucket_size) {
                min_value.push(page[0].fragment_mz);
                page.sort_by_key(|f| f.peptide_index);
            }
            let db = IndexedDatabase {
                peptides: peptides.clone(),
                fragment_mzs: ordered.iter().map(|f| f.fragment_mz).collect(),
                fragment_peptide_indices: ordered.iter().map(|f| f.peptide_index).collect(),
                min_value,
                bucket_size,
                ..Default::default()
            };
            for (mass, tol) in [
                (500.0, Tolerance::Da(0.0, 0.0)),
                (600.0, Tolerance::Da(0.0, 0.0)),
                (550.0, Tolerance::Da(-50.0, 50.0)),
                (550.0, Tolerance::Da(-1.0, 1.0)),
                (900.0, Tolerance::Da(0.0, 0.0)),
                (100.0, Tolerance::Da(-1.0, 1.0)),
                (1000.0, Tolerance::Da(-1.0, 1.0)),
                (700.0, Tolerance::Da(-500.0, 500.0)),
                (600.0, Tolerance::Ppm(-10.0, 10.0)),
            ] {
                let (lo, hi) = tol.bounds(mass);
                let mut expected = Vec::new();
                for &(center, tolerance) in &windows {
                    let (frag_lo, frag_hi) = tolerance.bounds(center);
                    for f in &fragments {
                        let precursor = db[f.peptide_index].monoisotopic;
                        if precursor >= lo
                            && precursor <= hi
                            && f.fragment_mz >= frag_lo
                            && f.fragment_mz <= frag_hi
                        {
                            expected.push((f.peptide_index, f.fragment_mz.to_bits()));
                        }
                    }
                }
                let query = db.query(mass, tol);
                let mut actual = Vec::new();
                query.page_search_batch(&windows, |f| {
                    actual.push((f.peptide_index, f.fragment_mz.to_bits()))
                });
                let mut unbatched: Vec<_> = windows
                    .iter()
                    .flat_map(|&(mass, tol)| query.page_search(mass, tol))
                    .map(|f| (f.peptide_index, f.fragment_mz.to_bits()))
                    .collect();
                expected.sort_unstable();
                actual.sort_unstable();
                unbatched.sort_unstable();
                assert_eq!(
                    actual, expected,
                    "bucket={bucket_size}, mass={mass}, tol={tol:?}"
                );
                assert_eq!(unbatched, expected);
            }
            db.query(600.0, Tolerance::Da(-1000.0, 1000.0))
                .page_search_batch(&[], |_| panic!("empty window list matched"));
        }
        IndexedDatabase::default()
            .query(600.0, Tolerance::Da(-1000.0, 1000.0))
            .page_search_batch(&windows, |_| panic!("empty database matched"));
    }

    #[cfg(target_arch = "x86_64")]
    #[test]
    fn avx2_scan_preserves_scalar_matches_order_and_tails() {
        if !std::arch::is_x86_feature_detected!("avx2") {
            return;
        }
        let values: Vec<_> = (0..80)
            .map(|i| match i % 11 {
                0 => f32::NAN,
                1 => f32::NEG_INFINITY,
                2 => f32::INFINITY,
                3 => -0.0,
                4 => 0.0,
                _ => (i % 11) as f32,
            })
            .collect();
        for offset in 0..8 {
            for len in 0..=65 {
                let masses = &values[offset..offset + len];
                for (lo, hi) in [
                    (-0.0, 0.0),
                    (5.0, 8.0),
                    (-f32::INFINITY, f32::INFINITY),
                    (f32::NAN, 8.0),
                    (0.0, f32::NAN),
                    (9.0, 2.0),
                ] {
                    let mut expected = Vec::new();
                    scan_masses_scalar(masses, lo, hi, |i| expected.push(i));
                    let mut actual = Vec::new();
                    // SAFETY: checked AVX2 support above.
                    unsafe {
                        scan_masses_avx2(masses, lo, hi, |i| actual.push(i));
                    }
                    assert_eq!(actual, expected, "offset={offset}, len={len}");
                }
            }
        }
    }

    #[test]
    fn binary_search_slice_smoke() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0];
        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.75, 3.5);
        assert_eq!(bounds, (1, 6));
        assert!(data[bounds.0] <= 1.75);
        assert_eq!(&data[bounds.0..bounds.1], &[1.5, 2.0, 2.5, 3.0, 3.5]);

        let bounds = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 0.0, 5.0);
        assert_eq!(bounds, (0, data.len()));
    }

    #[test]
    fn binary_search_slice_run() {
        // Make sure that our query returns the maximal set of indices
        let data = [1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0, 3.5, 4.0];
        let (left, right) = binary_search_slice(&data, |a: &f64, b| a.total_cmp(b), 1.5, 3.25);
        assert!(data[left] <= 1.5);
        assert!(data[right] > 3.25);
        assert_eq!(
            &data[left..right],
            &[1.0, 1.5, 1.5, 1.5, 1.5, 2.0, 2.5, 3.0, 3.0]
        );
    }

    #[test]
    fn digestion() {
        let fasta = r#"
        >sp|AAAAA
        MEWKLEQSMREQALLKAQLTQLK
        >sp|BBBBB
        RMEWKLEQSMREQALLKAQLTQLK
        "#;

        let fasta = Fasta::parse(fasta.into(), "rev_", false);

        // Make sure that FASTA parsed OK
        assert_eq!(
            fasta.targets,
            vec![
                (
                    Arc::from("sp|AAAAA".to_string()),
                    "MEWKLEQSMREQALLKAQLTQLK".into()
                ),
                (
                    Arc::from("sp|BBBBB".to_string()),
                    "RMEWKLEQSMREQALLKAQLTQLK".into()
                ),
            ]
        );

        let params = Parameters {
            bucket_size: 128,
            enzyme: EnzymeBuilder {
                missed_cleavages: Some(1),
                min_len: Some(6),
                max_len: Some(10),
                ..Default::default()
            },
            peptide_min_mass: 150.0,
            peptide_max_mass: 5000.0,
            ion_kinds: vec![Kind::B, Kind::Y],
            min_ion_index: 2,
            static_mods: HashMap::default(),
            variable_mods: [(ModificationSpecificity::ProteinN(None), vec![42.0])]
                .into_iter()
                .collect(),
            max_variable_mods: 2,
            decoy_tag: "rev_".into(),
            generate_decoys: false,
            fasta: "none".into(),
            prefilter: false,
            prefilter_chunk_size: 0,
            prefilter_low_memory: true,
        };

        let peptides = params.digest(&fasta);

        let expected = [
            "EQALLK",
            "LEQSMR",
            "AQLTQLK",
            "MEWKLEQSMR",
            "[+42]-MEWKLEQSMR",
        ]
        .into_iter()
        .map(String::from)
        .collect::<Vec<_>>();

        let sequences = peptides.iter().map(|p| p.to_string()).collect::<Vec<_>>();
        assert_eq!(expected, sequences);

        // All peptides are shared except for the protein N-term mod
        for peptide in &peptides[..4] {
            assert_eq!(peptide.proteins.len(), 2, "{:?}", peptide);
        }
        // Ensure that this mod is uniquely called as the first protein
        assert_eq!(
            peptides.last().unwrap().proteins,
            vec!["sp|AAAAA".to_string().into()]
        );
    }

    /// The interleaved resolver must agree exactly with the straightforward
    /// per-page `partition_point` pair it replaced, including at page
    /// boundaries, on partial final pages, and when the requested peptide-ID
    /// range is empty or covers the whole page.
    #[test]
    fn batched_page_bounds_match_partition_point() {
        fn reference(
            all_ids: &[PeptideIx],
            bucket_size: usize,
            total: usize,
            page: usize,
            lo: usize,
            hi: usize,
        ) -> (u32, u32) {
            let left = page * bucket_size;
            let right = ((page + 1) * bucket_size).min(total);
            let ids = &all_ids[left..right];
            let inner_left = ids.partition_point(|id| (id.0 as usize) < lo);
            let inner_right =
                inner_left + ids[inner_left..].partition_point(|id| (id.0 as usize) < hi);
            (inner_left as u32, inner_right as u32)
        }

        // Sorted-per-page ID columns with duplicates and gaps, plus a
        // deliberately short final page.
        for &bucket_size in &[1usize, 2, 3, 8, 17] {
            for &n_pages in &[1usize, 2, 5, 9] {
                for &tail in &[0usize, 1] {
                    let total = n_pages * bucket_size + tail;
                    if total == 0 {
                        continue;
                    }
                    let mut all_ids = Vec::with_capacity(total);
                    for page in 0..=n_pages {
                        let start = page * bucket_size;
                        let end = (start + bucket_size).min(total);
                        for k in start..end {
                            // Ascending within each page, with duplicates.
                            all_ids.push(PeptideIx(((k - start) as u32 / 2) * 3));
                        }
                    }
                    let pages = total.div_ceil(bucket_size);
                    let max_id = 3 * (bucket_size as u32 / 2 + 2);

                    for lo in 0..=max_id as usize {
                        for hi in lo..=max_id as usize {
                            let groups: Vec<(usize, usize, usize)> =
                                (0..pages).map(|p| (p, 0, 0)).collect();
                            let mut out = Vec::new();
                            resolve_page_peptide_bounds(
                                &all_ids, bucket_size, total, &groups, lo, hi, &mut out,
                            );
                            assert_eq!(out.len(), pages);
                            for (idx, &(page, _, _)) in groups.iter().enumerate() {
                                assert_eq!(
                                    out[idx],
                                    reference(&all_ids, bucket_size, total, page, lo, hi),
                                    "bucket_size={} total={} page={} lo={} hi={}",
                                    bucket_size, total, page, lo, hi
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}
