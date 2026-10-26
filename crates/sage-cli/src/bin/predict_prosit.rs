//! Offline Prosit 2023 timsTOF inference directly over a stored compact table.
use anyhow::{ensure, Context, Result};
use clap::{Arg, Command};
use sage_cli::native_prosit::{self, engine::TensorRt, Options};
use sage_core::{
    database::Builder,
    intensity_codec::IntensityEncoding,
    peptide_view::PeptideView,
    prosit::{ModificationPolicy, Tokenizer},
};
use sha2::{Digest, Sha256};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;
use std::time::Instant;

fn hash(path: &Path) -> Result<String> {
    let mut file = File::open(path)?;
    let mut hash = Sha256::new();
    let mut buffer = [0; 65536];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        hash.update(&buffer[..n]);
    }
    Ok(format!("{:x}", hash.finalize()))
}

fn write_json(path: &Path, value: &serde_json::Value) -> Result<()> {
    let mut file = File::create(path)?;
    serde_json::to_writer_pretty(&mut file, value)?;
    file.write_all(b"\n")?;
    file.sync_all()?;
    Ok(())
}

fn main() -> Result<()> {
    let wall = Instant::now();
    let matches = Command::new("predict_prosit")
        .about("Stream compact peptide rows through native Prosit 2023 TensorRT inference into Sage-compatible mmappet arrays")
        .arg(Arg::new("fasta").long("fasta").required(true))
        .arg(Arg::new("config").long("config").required(true).help("Same database JSON used to create the compact table"))
        .arg(Arg::new("compact").long("compact-peptides").required(true))
        .arg(Arg::new("engine").long("engine").required(true).help("Prosit 2023 timsTOF TensorRT 8.6 plan: Dense FP16/FP32 [N,29,6] or legacy FP16 [N,174]"))
        .arg(Arg::new("library").long("library").required(true).help("libprosit_trt.so built against the matching TensorRT runtime"))
        .arg(Arg::new("output").long("output").required(true).help("New directory; published only after complete, durable output"))
        .arg(Arg::new("min-charge").long("min-charge").default_value("2").value_parser(clap::value_parser!(usize)))
        .arg(Arg::new("max-charge").long("max-charge").default_value("4").value_parser(clap::value_parser!(usize)))
        .arg(Arg::new("ce").long("collision-energy").default_value("30").value_parser(clap::value_parser!(f32)))
        .arg(Arg::new("batch").long("batch-size").default_value("1000").value_parser(clap::value_parser!(usize)))
        .arg(Arg::new("top-k").long("top-k").default_value("12").value_parser(clap::value_parser!(usize)).help("Keep at most K most intense valid fragments per peptide/charge; 0 keeps all. Ties prefer earlier model fragment indices."))
        .arg(Arg::new("encoding").long("intensity-encoding").default_value("fp16").value_parser(["fp16", "sqrt8_basepeak_v1"]).help("FP16 preserves amplitudes; sqrt8_basepeak_v1 stores lossy normalized UINT8 values without scales"))
        .arg(Arg::new("start").long("row-start").default_value("0").value_parser(clap::value_parser!(usize)))
        .arg(Arg::new("count").long("row-count").value_parser(clap::value_parser!(usize)).help("Testing subset; partial exports intentionally omit Sage's full-table fingerprint"))
        .arg(Arg::new("modifications").long("modifications").default_value("original").value_parser(["original", "all"]))
        .arg(Arg::new("unsupported").long("unsupported").default_value("error").value_parser(["error", "skip"]))
        .arg(Arg::new("unimod").long("unimod-db-path"))
        .arg(Arg::new("reference").long("reference-tsv").help("Optional validation only: materialize sequences to a new TSV outside the inference timer"))
        .get_matches();
    let string = |name: &str| matches.get_one::<String>(name).unwrap().as_str();
    let number = |name: &str| *matches.get_one::<usize>(name).unwrap();
    if let Some(path) = matches.get_one::<String>("unimod") {
        sage_core::unimod::set_active_table_from_path(Path::new(path))
            .map_err(anyhow::Error::msg)?;
    }
    let output = Path::new(string("output"));
    ensure!(
        !output.exists(),
        "output already exists: {}",
        output.display()
    );
    let parent = output
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent)?;
    let name = output
        .file_name()
        .context("output must name a directory")?
        .to_string_lossy();
    let staging = parent.join(format!(".{name}.partial-{}", std::process::id()));
    ensure!(
        !staging.exists(),
        "staging directory exists: {}",
        staging.display()
    );

    let mut config: serde_json::Value = serde_json::from_slice(&std::fs::read(string("config"))?)?;
    sage_cli::input::resolve_unimod_refs_in_database(&mut config)?;
    let mut builder: Builder = serde_json::from_value(config)?;
    builder.fasta = Some(string("fasta").to_owned());
    let parameters = builder.make_parameters();
    let fasta = sage_cloudpath::util::read_fasta(
        &parameters.fasta,
        &parameters.decoy_tag,
        parameters.generate_decoys,
    )?;
    let table =
        sage_cloudpath::compact_store::open(Path::new(string("compact")), &parameters, &fasta)?;
    let start = number("start");
    ensure!(start <= table.rows.len(), "row-start outside table");
    let options = Options {
        row_start: start,
        row_count: matches
            .get_one::<usize>("count")
            .copied()
            .unwrap_or(table.rows.len() - start),
        min_charge: number("min-charge"),
        max_charge: number("max-charge"),
        collision_energy: *matches.get_one::<f32>("ce").unwrap(),
        batch_size: number("batch"),
        top_k: number("top-k"),
        intensity_encoding: if string("encoding") == "fp16" {
            IntensityEncoding::Fp16
        } else {
            IntensityEncoding::Sqrt8Basepeak
        },
        skip_unsupported: string("unsupported") == "skip",
    };
    options.validate(&table)?;
    let policy = if string("modifications") == "original" {
        ModificationPolicy::Original
    } else {
        ModificationPolicy::StripAll
    };
    let mut modifications = HashMap::new();
    for mass in parameters
        .static_mods
        .values()
        .chain(parameters.variable_mods.values().flatten())
    {
        if let Some(id) = sage_core::unimod::lookup_reverse(*mass) {
            modifications.insert(mass.to_bits(), id);
        }
    }
    let tokenizer = Tokenizer::new(policy, modifications);
    if let Some(path) = matches.get_one::<String>("reference") {
        let file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path)?;
        let mut writer = BufWriter::new(file);
        writeln!(
            writer,
            "peptide_row\tlength\tdecoy\tsupported\tsequence\ttokens"
        )?;
        for row in start..start + options.row_count {
            let view = table.view(row);
            let mut tokens = [0; 30];
            let supported = tokenizer.encode(&view, &mut tokens).is_ok();
            let peptide = table.materialize(&table.rows.get(row));
            let token_text = tokens
                .iter()
                .map(i32::to_string)
                .collect::<Vec<_>>()
                .join(",");
            writeln!(
                writer,
                "{row}\t{}\t{}\t{supported}\t{peptide}\t{token_text}",
                view.len(),
                view.decoy()
            )?;
        }
        writer.flush()?;
    }
    let engine_path = Path::new(string("engine"));
    let library_path = Path::new(string("library"));
    let engine_sha256 = hash(engine_path)?;
    let library_sha256 = hash(library_path)?;
    let compact_meta_sha256 = hash(&Path::new(string("compact")).join("meta.json"))?;
    let mut engine = TensorRt::open(library_path, engine_path, options.batch_size)?;
    let setup_seconds = wall.elapsed().as_secs_f64();
    eprintln!(
        "Predicting {} peptide rows × {} charges; output {}",
        options.row_count,
        options.max_charge - options.min_charge + 1,
        output.display()
    );
    let metrics = native_prosit::run(&table, &tokenizer, &mut engine, &options, &staging)
        .with_context(|| {
            format!(
                "incomplete files, if any, remain unpublished at {}",
                staging.display()
            )
        })?;
    ensure!(
        metrics.fragments > 0,
        "no fragments retained; output remains unpublished at {}",
        staging.display()
    );
    let complete_table = start == 0 && options.row_count == table.rows.len();
    let mut metadata = serde_json::json!({
        "format": "sage_native_prosit", "version": if options.intensity_encoding == IntensityEncoding::Fp16 { 1 } else { 2 },
        "intensity_encoding": options.intensity_encoding.as_str(),
        "normalization": if options.intensity_encoding == IntensityEncoding::Fp16 { "source_amplitude" } else { "per_spectrum_base_peak" },
        "rounding": if options.intensity_encoding == IntensityEncoding::Fp16 { "source_fp16" } else { "nearest_even" },
        "model": "Prosit_2023_intensity_timsTOF", "engine_sha256": engine_sha256,
        "native_bridge_sha256": library_sha256, "compact_meta_sha256": compact_meta_sha256,
        "source_fingerprint": table.stored_fingerprint,
        "row_start": start, "row_count": options.row_count, "total_peptide_rows": table.rows.len(),
        "min_charge": options.min_charge, "max_charge": options.max_charge,
        "collision_energy": options.collision_energy, "modifications": string("modifications"),
        "unsupported": string("unsupported"), "batch_size": options.batch_size,
        "top_k": options.top_k,
        "top_k_selection": if engine.dense_output() { "after validity filtering; descending Dense intensity, ties by ascending model fragment index; retained output in intensity order" } else { "after validity filtering; descending FP16 intensity, ties by ascending model fragment index; retained output in model order" },
        "model_output": engine.model_output(),
        "fragment_order": if engine.dense_output() { "descending_intensity" } else { "model" },
        "slot_order": "(peptide_row - row_start) * charge_span + charge - min_charge",
        "annotation_encoding": "sage_backbone_b87_y87_v1",
        "postprocessing": "gpu", "native_bridge_abi": 3,
        "filter": "ordinal < peptide_length; fragment_charge <= precursor_charge; intensity >= 0; retain zero",
        "metrics": metrics, "setup_seconds": setup_seconds,
        "wall_seconds_before_publication": wall.elapsed().as_secs_f64(),
        "reuse_existing_predictions": false
    });
    if complete_table {
        metadata["dumped_peptides_sha256"] =
            serde_json::json!(table
                .stored_fingerprint
                .as_ref()
                .context("compact table lacks fingerprint")?);
    }
    write_json(&staging.join("index.mmappet/meta.json"), &metadata)?;
    write_json(&staging.join("arrays.mmappet/meta.json"), &metadata)?;
    write_json(&staging.join("manifest.json"), &metadata)?;
    for directory in [
        staging.join("arrays.mmappet"),
        staging.join("index.mmappet"),
        staging.clone(),
    ] {
        File::open(directory)?.sync_all()?;
    }
    // A single producer owns this new output. Refuse to replace an existing dataset.
    ensure!(
        !output.exists(),
        "output appeared during inference; refusing publication"
    );
    std::fs::rename(&staging, output)?;
    File::open(parent)?.sync_all()?;
    println!("{}", serde_json::to_string(&metadata)?);
    Ok(())
}
