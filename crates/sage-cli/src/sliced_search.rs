//! Searching spectra against a `CompactPeptides` table one precursor-mass slice at
//! a time, so only one slice's peptides and fragment index are ever in memory. See
//! `docs/ai/sliced_search.md`.

use crate::process_memory::available_memory_bytes;
use crate::runner::{Predictions, Runner};
use log::info;
use rayon::prelude::*;
use sage_core::compact::{drop_in_parallel, CompactPeptide, CompactPeptides};
use sage_core::database::{IndexedDatabase, PeptideIx};
use sage_core::scoring::{Feature, RankingScore, Scorer};
use sage_core::slicing::SlicePlan;
use sage_core::spectrum::ProcessedSpectrum;
use std::time::Instant;

/// While a slice's index is built, a row costs its `Peptide` plus, per fragment,
/// the index entry and the copies sorting and unzipping make.
const BYTES_PER_ROW: usize = 200;
const BYTES_PER_FRAGMENT: usize = 16;
const AVAILABLE_MEMORY_PER_SLICE: f64 = 0.5;
/// Added to the widest candidate window, against float rounding at slice edges.
const WINDOW_MARGIN_DA: f32 = 0.01;
/// Caps each slice's fragments instead of sizing slices by available memory.
const SLICE_FRAGMENTS_VARIABLE: &str = "SAGE_SLICE_FRAGMENTS";

/// One spectrum searched at one precursor charge.
struct Query {
    spectrum: usize,
    charge: u8,
    /// Lowest candidate mass; decides the slice.
    window_low: f32,
    /// The spectrum's own precursor isn't annotated with `charge`, so it is scored
    /// through a copy that is.
    needs_annotation: bool,
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

/// Search `spectra` slice by slice; returns the same features, in the same order,
/// as one search over the whole table, with `peptide_idx` set to table rows.
pub fn search(
    runner: &Runner,
    table: &CompactPeptides,
    predictions: &Predictions,
    spectra: &[ProcessedSpectrum],
) -> anyhow::Result<Vec<Feature>> {
    let start = Instant::now();
    let no_peptides = IndexedDatabase::default();
    let planner = runner.make_scorer(&no_peptides, 0, predictions);
    let queries = plan_queries(runner, &planner, spectra);
    let widest_window = queries
        .iter()
        .map(|query| {
            let precursor = &spectra[query.spectrum].precursors[0];
            let (low, high) = planner.candidate_mass_window(precursor, query.charge);
            high - low
        })
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
        "sliced search: {} queries over {} slices (widest window {:.3} Da)",
        queries.len(),
        plan.slices.len(),
        widest_window
    );

    let mut features_by_query: Vec<Vec<Feature>> = Vec::new();
    features_by_query.resize_with(queries.len(), Vec::new);
    for (slice_index, (slice, query_indices)) in
        plan.slices.iter().zip(&queries_by_slice).enumerate()
    {
        if query_indices.is_empty() {
            continue;
        }
        let build_start = Instant::now();
        let peptides = table.materialize_rows(slice.rows.clone());
        let db = runner
            .parameters
            .database
            .clone()
            .build_from_peptides(peptides);
        let build_seconds = build_start.elapsed().as_secs_f64();

        let search_start = Instant::now();
        let scorer = Scorer {
            override_precursor_charge: false,
            ..runner.make_scorer(&db, slice.rows.start, predictions)
        };
        let first_row = slice.rows.start as u32;
        let scored: Vec<(usize, Vec<Feature>)> = query_indices
            .par_iter()
            .map(|&index| {
                let query = &queries[index];
                let mut features = score_query(&scorer, &spectra[query.spectrum], query);
                for feature in &mut features {
                    feature.peptide_idx = PeptideIx(feature.peptide_idx.0 + first_row);
                }
                (index, features)
            })
            .collect();
        info!(
            "- slice {}/{}: {} rows, {} fragments, {} queries; build {:.1} s, search {:.1} s",
            slice_index + 1,
            plan.slices.len(),
            slice.rows.len(),
            db.size(),
            query_indices.len(),
            build_seconds,
            search_start.elapsed().as_secs_f64(),
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
        "- sliced search: {:8} ms ({} spectra)",
        start.elapsed().as_millis(),
        spectra.len()
    );
    Ok(features)
}

/// One query per spectrum and precursor charge a whole-table search would try,
/// for the spectra it would search at all, in spectrum order.
fn plan_queries(runner: &Runner, planner: &Scorer, spectra: &[ProcessedSpectrum]) -> Vec<Query> {
    let mut queries = Vec::new();
    for (index, spectrum) in spectra.iter().enumerate() {
        if spectrum.peaks.len() < runner.parameters.min_peaks || spectrum.level != 2 {
            continue;
        }
        let precursor = spectrum.precursors.first().unwrap_or_else(|| {
            panic!("missing MS1 precursor for {}", spectrum.id);
        });
        for charge in planner.candidate_charges(precursor) {
            queries.push(Query {
                spectrum: index,
                charge,
                window_low: planner.candidate_mass_window(precursor, charge).0,
                needs_annotation: precursor.charge != Some(charge),
            });
        }
    }
    queries
}

fn score_query(scorer: &Scorer, spectrum: &ProcessedSpectrum, query: &Query) -> Vec<Feature> {
    if !query.needs_annotation {
        return scorer.score(spectrum);
    }
    let mut annotated = spectrum.clone();
    annotated.precursors[0].charge = Some(query.charge);
    scorer.score(&annotated)
}

/// Every spectrum's features, spectra in query order. A spectrum searched at
/// several charges gets its charges' features ranked together.
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
            .is_some_and(|next| next.spectrum == queries[index].spectrum);
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
