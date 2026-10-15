//! Searching spectra against a `CompactPeptides` table one precursor-mass slice at
//! a time, so only one slice's peptides and fragment index are ever in memory. Spectra
//! are read by index when their slice is searched, so they aren't held either. See
//! `docs/ai/sliced_search.md` and `docs/ai/streaming_spectra.md`.

use crate::indexed_spectra::{IndexedReader, IndexedSpectra};
use crate::process_cpu::{cpu_seconds, utilisation};
use crate::process_memory::available_memory_bytes;
use crate::runner::{Predictions, Runner};
use log::info;
use rayon::prelude::*;
use sage_core::compact::{drop_in_parallel, CompactPeptide, CompactPeptides};
use sage_core::database::{IndexedDatabase, PeptideIx};
use sage_core::scoring::{Feature, RankingScore, Scorer};
use sage_core::slicing::{rows_reaching, SlicePlan};
use sage_core::spectrum::{Precursor, ProcessedSpectrum, SpectrumProcessor};
use std::time::Instant;

/// While a slice's index is built, a row costs its `Peptide` plus, per fragment,
/// the index entry and the copies sorting and unzipping make.
const BYTES_PER_ROW: usize = 200;
const BYTES_PER_FRAGMENT: usize = 16;
const AVAILABLE_MEMORY_PER_SLICE: f64 = 0.5;
/// Added to the candidate windows, against float rounding at slice edges.
const WINDOW_MARGIN_DA: f32 = 0.01;
/// Caps each slice's fragments instead of sizing slices by available memory.
const SLICE_FRAGMENTS_VARIABLE: &str = "SAGE_SLICE_FRAGMENTS";

/// One spectrum searched at one precursor charge.
struct Query {
    file: usize,
    spectrum: usize,
    charge: u8,
    /// Lowest candidate mass; decides the slice.
    window_low: f32,
    window_high: f32,
    /// The spectrum as read would search other charges too (or not this one), so it is
    /// scored through a copy whose only precursor carries `charge`.
    needs_annotation: bool,
}

impl Query {
    fn same_spectrum(&self, other: &Query) -> bool {
        self.file == other.file && self.spectrum == other.spectrum
    }
}

/// How big one slice may be.
enum SliceBudget {
    Fragments(usize),
    Bytes(usize),
}

impl SliceBudget {
    fn from_environment() -> anyhow::Result<Self> {
        if let Ok(fragments) = std::env::var(SLICE_FRAGMENTS_VARIABLE) {
            return Ok(Self::Fragments(fragments.parse()?));
        }
        let available = available_memory_bytes()?;
        Ok(Self::Bytes(
            (available as f64 * AVAILABLE_MEMORY_PER_SLICE) as usize,
        ))
    }

    fn limit(&self) -> usize {
        match self {
            Self::Fragments(fragments) => *fragments,
            Self::Bytes(bytes) => *bytes,
        }
    }

    fn cost(&self, table: &CompactPeptides, row: &CompactPeptide) -> usize {
        match self {
            Self::Fragments(_) => table.fragment_count(row),
            Self::Bytes(_) => BYTES_PER_ROW + BYTES_PER_FRAGMENT * table.fragment_count(row),
        }
    }
}

/// Search the spectra of `files` slice by slice; returns the same features, in the
/// same order (files, then spectra), as one search over the whole table, with
/// `peptide_idx` set to table rows.
pub fn search(
    runner: &Runner,
    table: &CompactPeptides,
    predictions: &Predictions,
    files: &[IndexedSpectra],
    processor: &SpectrumProcessor,
) -> anyhow::Result<Vec<Feature>> {
    let start = Instant::now();
    let no_peptides = IndexedDatabase::default();
    let planner = runner.make_scorer(&no_peptides, 0, predictions);
    let queries = plan_queries(&planner, files);
    let widest_window = queries
        .iter()
        .map(|query| query.window_high - query.window_low)
        .fold(0.0_f32, f32::max);

    let budget = SliceBudget::from_environment()?;
    let plan = SlicePlan::new(
        &table.rows,
        |row| row.monoisotopic,
        |row| budget.cost(table, row),
        budget.limit(),
        widest_window + WINDOW_MARGIN_DA,
    );
    let mut queries_by_slice = vec![Vec::new(); plan.slices.len()];
    for (index, query) in queries.iter().enumerate() {
        queries_by_slice[plan.slice_for(query.window_low)].push(index);
    }
    info!(
        "sliced search: {} queries over {} slices (widest window {:.3} Da), planned in {:.1} s",
        queries.len(),
        plan.slices.len(),
        widest_window,
        start.elapsed().as_secs_f64()
    );

    let mut features_by_query: Vec<Vec<Feature>> = Vec::new();
    features_by_query.resize_with(queries.len(), Vec::new);
    for (slice_index, query_indices) in queries_by_slice.iter().enumerate() {
        if query_indices.is_empty() {
            continue;
        }
        let rows = rows_reached_by(table, &queries, query_indices);
        let build = Phase::start();
        let peptides = table.materialize_rows(rows.clone());
        let db = runner
            .parameters
            .database
            .clone()
            .build_from_peptides(peptides);
        let build = build.finish();

        let search = Phase::start();
        let scorer = Scorer {
            override_precursor_charge: false,
            ..runner.make_scorer(&db, rows.start, predictions)
        };
        let first_row = rows.start as u32;
        let same_spectrum_runs: Vec<&[usize]> = query_indices
            .chunk_by(|&a, &b| queries[a].same_spectrum(&queries[b]))
            .collect();
        let scored: Vec<(usize, Vec<Feature>)> = same_spectrum_runs
            .par_iter()
            .map_init(
                || files.iter().map(IndexedSpectra::reader).collect::<Vec<_>>(),
                |readers, run| {
                    score_spectrum(runner, &scorer, files, processor, readers, &queries, run)
                },
            )
            .flatten_iter()
            .map(|(index, mut features)| {
                for feature in &mut features {
                    feature.peptide_idx = PeptideIx(feature.peptide_idx.0 + first_row);
                }
                (index, features)
            })
            .collect();
        let search = search.finish();
        info!(
            "- slice {}/{}: {} rows, {} fragments, {} queries; build {:.1} s at {:.0}% CPU, search {:.1} s at {:.0}% CPU",
            slice_index + 1,
            plan.slices.len(),
            rows.len(),
            db.size(),
            query_indices.len(),
            build.wall_seconds,
            100.0 * build.utilisation,
            search.wall_seconds,
            100.0 * search.utilisation,
        );
        for (index, features) in scored {
            features_by_query[index] = features;
        }
        drop_in_parallel(db.peptides);
    }

    let features = in_spectrum_order(
        &queries,
        features_by_query,
        runner.parameters.report_psms,
        runner.parameters.ranking_score,
    );
    info!(
        "- sliced search: {:8} ms ({} queries)",
        start.elapsed().as_millis(),
        queries.len()
    );
    Ok(features)
}

/// Wall time and how busy it kept the cores.
struct Phase {
    wall_seconds: f64,
    utilisation: f64,
}

struct PhaseStart {
    wall: Instant,
    cpu_seconds: Option<f64>,
}

impl Phase {
    fn start() -> PhaseStart {
        PhaseStart {
            wall: Instant::now(),
            cpu_seconds: cpu_seconds(),
        }
    }
}

impl PhaseStart {
    fn finish(self) -> Phase {
        let wall_seconds = self.wall.elapsed().as_secs_f64();
        let cpu = match (self.cpu_seconds, cpu_seconds()) {
            (Some(before), Some(after)) => after - before,
            _ => 0.0,
        };
        Phase {
            wall_seconds,
            utilisation: utilisation(cpu, wall_seconds),
        }
    }
}

/// The rows a slice must hold for its queries: from just below their lowest window
/// start to their highest window end (see `slicing::rows_reaching`).
fn rows_reached_by(
    table: &CompactPeptides,
    queries: &[Query],
    query_indices: &[usize],
) -> std::ops::Range<usize> {
    let low = query_indices
        .iter()
        .map(|&index| queries[index].window_low)
        .fold(f32::INFINITY, f32::min);
    let high = query_indices
        .iter()
        .map(|&index| queries[index].window_high)
        .fold(f32::NEG_INFINITY, f32::max);
    rows_reaching(
        &table.rows,
        |row| row.monoisotopic,
        low - WINDOW_MARGIN_DA,
        high + WINDOW_MARGIN_DA,
    )
}

/// One query per spectrum and precursor charge a whole-table search would try, for the
/// MS2 spectra of `files`, in file then spectrum order. Read from headers only; the
/// `min_peaks` check waits for scoring (a spectrum below it yields no features either
/// way).
fn plan_queries(planner: &Scorer, files: &[IndexedSpectra]) -> Vec<Query> {
    let mut queries = Vec::new();
    for (file, spectra) in files.iter().enumerate() {
        let headers: Vec<_> = (0..spectra.len())
            .into_par_iter()
            .map_init(
                || spectra.reader(),
                |reader, index| spectra.header(index, reader),
            )
            .collect();
        for (index, header) in headers.into_iter().enumerate() {
            let Some(header) = header.filter(|header| header.level == 2) else {
                continue;
            };
            let Some(precursor) = header.precursors.first() else {
                panic!("missing MS1 precursor for spectrum {index} of file {file}");
            };
            let charges = planner.candidate_charges(&header.precursors);
            for charge in charges.iter() {
                let (window_low, window_high) = planner.candidate_mass_window(precursor, charge);
                queries.push(Query {
                    file,
                    spectrum: index,
                    charge,
                    window_low,
                    window_high,
                    needs_annotation: !charges.is_only(charge),
                });
            }
        }
    }
    queries
}

/// Reads the spectrum of a run of queries once and scores each of its charges.
fn score_spectrum<'a>(
    runner: &Runner,
    scorer: &Scorer,
    files: &'a [IndexedSpectra],
    processor: &SpectrumProcessor,
    readers: &mut [IndexedReader<'a>],
    queries: &[Query],
    run: &[usize],
) -> Vec<(usize, Vec<Feature>)> {
    let first = &queries[run[0]];
    let spectrum = files[first.file].processed(first.spectrum, processor, &mut readers[first.file]);
    run.iter()
        .map(|&index| {
            let features = match &spectrum {
                Some(spectrum) if runner.is_searched(spectrum) => {
                    score_query(scorer, spectrum, &queries[index])
                }
                _ => Vec::new(),
            };
            (index, features)
        })
        .collect()
}

fn score_query(scorer: &Scorer, spectrum: &ProcessedSpectrum, query: &Query) -> Vec<Feature> {
    if !query.needs_annotation {
        return scorer.score(spectrum);
    }
    let mut annotated = spectrum.clone();
    annotated.precursors = vec![Precursor {
        charge: Some(query.charge),
        ..spectrum.precursors[0].clone()
    }];
    scorer.score(&annotated)
}

/// Every spectrum's features, spectra in query order. A spectrum searched at several
/// charges gets its charges' features ranked together.
fn in_spectrum_order(
    queries: &[Query],
    features_by_query: Vec<Vec<Feature>>,
    report_psms: usize,
    ranking_score: RankingScore,
) -> Vec<Feature> {
    let mut features = Vec::new();
    let mut same_spectrum: Vec<Vec<Feature>> = Vec::new();
    for (index, query_features) in features_by_query.into_iter().enumerate() {
        same_spectrum.push(query_features);
        let spectrum_continues = queries
            .get(index + 1)
            .is_some_and(|next| next.same_spectrum(&queries[index]));
        if spectrum_continues {
            continue;
        }
        match same_spectrum.len() {
            1 => features.extend(same_spectrum.pop().unwrap()),
            _ => features.extend(rank_across_charges(
                std::mem::take(&mut same_spectrum),
                report_psms,
                ranking_score,
            )),
        }
    }
    features
}

/// One spectrum's PSMs from searches at several charges, ranked together: the best
/// `report_psms` by rank key, with `rank`, `delta_best` and `delta_next` recomputed
/// across charges. `delta_next` reaches past the reported PSMs through each charge's
/// best unreported candidate, which its last reported PSM's own `delta_next` gives
/// away.
fn rank_across_charges(
    by_charge: Vec<Vec<Feature>>,
    report_psms: usize,
    ranking_score: RankingScore,
) -> Vec<Feature> {
    let best_unreported = by_charge
        .iter()
        .filter_map(|features| features.last())
        .map(|last| last.rank_key(ranking_score) - last.delta_next)
        .fold(0.0_f64, f64::max);
    let mut merged: Vec<Feature> = by_charge.into_iter().flatten().collect();
    merged.sort_by(|a, b| {
        b.rank_key(ranking_score)
            .total_cmp(&a.rank_key(ranking_score))
    });
    let keys: Vec<f64> = merged.iter().map(|f| f.rank_key(ranking_score)).collect();
    merged.truncate(report_psms);
    for (index, feature) in merged.iter_mut().enumerate() {
        let next_reported = keys.get(index + 1).copied().unwrap_or_default();
        feature.rank = index as u32 + 1;
        feature.delta_best = keys[0] - keys[index];
        feature.delta_next = keys[index] - next_reported.max(best_unreported);
    }
    merged
}

/// Renumber `features` from table rows to a database holding only their peptides,
/// in table-row order, and return that database.
pub fn reported_database(table: &CompactPeptides, features: &mut [Feature]) -> IndexedDatabase {
    let mut rows: Vec<usize> = features.iter().map(|f| f.peptide_idx.0 as usize).collect();
    rows.par_sort_unstable();
    rows.dedup();
    features.par_iter_mut().for_each(|feature| {
        let row = feature.peptide_idx.0 as usize;
        let position = rows
            .binary_search(&row)
            .expect("every feature's row was collected");
        feature.peptide_idx = PeptideIx(position as u32);
    });
    info!("reporting {} distinct peptides", rows.len());
    table.database_for_rows(&rows)
}
