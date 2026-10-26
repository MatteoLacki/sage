//! Bounded preparation → GPU inference → native mmappet writing pipeline.
pub mod engine;
mod storage;

use anyhow::{anyhow, ensure, Result};
use sage_core::compact::CompactPeptides;
use sage_core::intensity_codec::IntensityEncoding;
use sage_core::prosit::{Tokenizer, SEQUENCE_LENGTH};
use serde::Serialize;
use std::collections::BTreeMap;
use std::path::Path;
use std::sync::mpsc::sync_channel;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct Options {
    pub row_start: usize,
    pub row_count: usize,
    pub min_charge: usize,
    pub max_charge: usize,
    pub collision_energy: f32,
    pub batch_size: usize,
    /// Maximum retained fragments per peptide/charge; zero retains all.
    pub top_k: usize,
    pub intensity_encoding: IntensityEncoding,
    pub skip_unsupported: bool,
}

impl Options {
    pub fn validate(&self, table: &CompactPeptides) -> Result<()> {
        ensure!(
            self.min_charge >= 1 && self.max_charge <= 6 && self.min_charge <= self.max_charge,
            "charges must be an ordered range within 1..6"
        );
        ensure!(
            self.collision_energy.is_finite() && self.collision_energy >= 0.0,
            "collision energy must be finite and nonnegative"
        );
        ensure!(
            (1..=1000).contains(&self.batch_size),
            "batch size must be 1..1000"
        );
        ensure!(
            self.row_start <= table.rows.len()
                && self.row_count <= table.rows.len() - self.row_start,
            "row range is outside the compact table"
        );
        self.slots()?;
        Ok(())
    }
    fn charge_span(&self) -> usize {
        self.max_charge - self.min_charge + 1
    }
    fn slots(&self) -> Result<usize> {
        self.row_count
            .checked_mul(self.charge_span())
            .ok_or_else(|| anyhow!("prediction slot count overflow"))
    }
}

struct Input {
    slot: usize,
}
struct Batch {
    first_slot: usize,
    slots: usize,
    tokens: Vec<i32>,
    charges: Vec<f32>,
    energies: Vec<f32>,
    output: engine::Predictions,
    inputs: Vec<Input>,
}
impl Batch {
    fn new(capacity: usize, top_k: usize, encoding: IntensityEncoding) -> Self {
        Self {
            first_slot: 0,
            slots: 0,
            tokens: Vec::with_capacity(capacity * 30),
            charges: Vec::with_capacity(capacity * 6),
            energies: Vec::with_capacity(capacity),
            output: engine::Predictions::new(top_k, encoding),
            inputs: Vec::with_capacity(capacity),
        }
    }
    fn clear(&mut self) {
        self.tokens.clear();
        self.charges.clear();
        self.energies.clear();
        self.inputs.clear();
        self.output.data.clear();
    }
}

#[derive(Default, Debug, Serialize)]
pub struct Metrics {
    pub peptide_rows: usize,
    pub prediction_slots: usize,
    pub predicted_slots: usize,
    pub skipped_peptide_rows: BTreeMap<String, usize>,
    pub fragments: i64,
    pub tokenize_seconds: f64,
    pub inference_seconds: f64,
    pub pack_write_seconds: f64,
    pub flush_seconds: f64,
    pub pipeline_seconds: f64,
}

pub fn run(
    table: &CompactPeptides,
    tokenizer: &Tokenizer,
    engine: &mut impl engine::Engine,
    options: &Options,
    staging: &Path,
) -> Result<Metrics> {
    options.validate(table)?;
    let writer = storage::Writer::new(staging, options.intensity_encoding)?;
    let total = options.slots()?;
    let (free_tx, free_rx) = sync_channel::<Batch>(3);
    let (prepared_tx, prepared_rx) = sync_channel::<Batch>(2);
    let (write_tx, write_rx) = sync_channel::<Batch>(2);
    for _ in 0..3 {
        free_tx.send(Batch::new(
            options.batch_size,
            options.top_k,
            options.intensity_encoding,
        ))?;
    }
    let started = Instant::now();
    let result = std::thread::scope(|scope| -> Result<Metrics> {
        let prepare = scope.spawn(move || -> Result<Metrics> {
            let mut metrics = Metrics { peptide_rows: options.row_count, prediction_slots: total, ..Metrics::default() };
            let mut last_row = None;
            let mut row_tokens = [0; SEQUENCE_LENGTH];
            let mut row_error = None;
            for start in (0..total).step_by(options.batch_size) {
                let mut batch = free_rx.recv().map_err(|_| anyhow!("output writer stopped"))?;
                let tick = Instant::now();
                batch.clear();
                batch.first_slot = start;
                batch.slots = options.batch_size.min(total - start);
                for slot in 0..batch.slots {
                    let row = options.row_start + (start + slot) / options.charge_span();
                    let charge = options.min_charge + (start + slot) % options.charge_span();
                    let view = table.view(row);
                    if last_row != Some(row) {
                        row_error = tokenizer.encode(&view, &mut row_tokens).err();
                        last_row = Some(row);
                        if let Some(error) = row_error {
                            ensure!(options.skip_unsupported, "unsupported peptide row {row}: {} (use --unsupported skip to write missing-prediction slots)", error.reason());
                            *metrics.skipped_peptide_rows.entry(error.reason().to_owned()).or_default() += 1;
                        }
                    }
                    if row_error.is_some() { continue; }
                    batch.tokens.extend_from_slice(&row_tokens);
                    let offset = batch.charges.len();
                    batch.charges.resize(offset + 6, 0.0);
                    batch.charges[offset + charge - 1] = 1.0;
                    batch.energies.push(options.collision_energy * 0.01_f32);
                    batch.inputs.push(Input { slot });
                }
                batch.output.resize(batch.inputs.len());
                metrics.predicted_slots += batch.inputs.len();
                metrics.tokenize_seconds += tick.elapsed().as_secs_f64();
                prepared_tx.send(batch).map_err(|_| anyhow!("inference stopped"))?;
            }
            Ok(metrics)
        });
        let write = scope.spawn(move || -> Result<(usize, i64, f64, f64)> {
            let mut writer = writer;
            let mut seconds = 0.0;
            for batch in write_rx {
                let tick = Instant::now();
                writer.append(&batch)?;
                seconds += tick.elapsed().as_secs_f64();
                // The producer may have finished already. No further reuse is needed then.
                let _ = free_tx.send(batch);
            }
            let tick = Instant::now();
            let (slots, fragments) = writer.finish()?;
            Ok((slots, fragments, seconds, tick.elapsed().as_secs_f64()))
        });
        let inference = (|| -> Result<f64> {
            let mut seconds = 0.0;
            for mut batch in &prepared_rx {
                if !batch.inputs.is_empty() {
                    let tick = Instant::now();
                    engine.predict(
                        &batch.tokens,
                        &batch.charges,
                        &batch.energies,
                        &mut batch.output,
                    )?;
                    seconds += tick.elapsed().as_secs_f64();
                }
                write_tx
                    .send(batch)
                    .map_err(|_| anyhow!("output writer stopped"))?;
            }
            Ok(seconds)
        })();
        drop(prepared_rx);
        drop(write_tx);
        let prepared = prepare
            .join()
            .map_err(|_| anyhow!("preparation thread panicked"))?;
        let written = write
            .join()
            .map_err(|_| anyhow!("writer thread panicked"))?;
        // Prefer the original validation/I/O failure over its resulting channel error.
        let (slots, fragments, writing, flushing) = written?;
        let inference_seconds = inference?;
        let mut metrics = prepared?;
        metrics.inference_seconds = inference_seconds;
        ensure!(slots == total, "incomplete positional index");
        metrics.fragments = fragments;
        metrics.pack_write_seconds = writing;
        metrics.flush_seconds = flushing;
        metrics.pipeline_seconds = started.elapsed().as_secs_f64();
        Ok(metrics)
    });
    result
}

#[cfg(test)]
mod tests;
