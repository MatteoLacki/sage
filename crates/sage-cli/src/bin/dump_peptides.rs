use clap::{Arg, Command, ValueHint};
use sage_core::database::Builder;

/// Digest a FASTA database and dump every target/decoy peptide, in the same
/// `Peptide::Display` inline-mod string format `sage` itself uses
/// internally, without building the fragment index or requiring any
/// spectra input. Used to generate the sequence list that RT/IIM
/// predictions (see `git/featureprediction`) are computed for.
///
/// Only depends on `--fasta` plus the same digestion-only settings `sage`
/// itself reads from its `database` config section (enzyme, static/variable
/// mods, mass bounds, decoy_tag) — unrelated search settings like
/// `precursor_tol`/`fragment_tol` play no part in digestion and are not
/// required here. Callers (e.g. the necroflow pipeline rule) are
/// responsible for slicing their full sage config down to just the
/// `database` subdictionary before passing it via `--config` — this binary
/// itself never reads or requires the rest of that config.
fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("dump_peptides")
        .version(clap::crate_version!())
        .about("Dump digested target/decoy peptides for a FASTA database, without searching spectra")
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .required(true)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to FASTA database")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("config")
                .short('c')
                .long("config")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Optional path to a JSON `database` section (enzyme, static/variable \
                     mods, mass bounds, decoy_tag — same shape as `sage`'s own config \
                     `database` key) to override defaults. Omit to use sage's defaults.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output")
                .short('o')
                .long("output")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .default_value("peptides.parquet")
                .help("Output parquet path")
                .value_hint(ValueHint::FilePath),
        )
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let output = matches.get_one::<String>("output").unwrap();

    let mut builder: Builder = match matches.get_one::<String>("config") {
        Some(path) => serde_json::from_slice(&std::fs::read(path)?)?,
        None => Builder::default(),
    };
    builder.fasta = Some(fasta_path.clone());
    let params = builder.make_parameters();

    let fasta = sage_cloudpath::util::read_fasta(
        &params.fasta,
        &params.decoy_tag,
        params.generate_decoys,
    )?;

    let decoy_tag = params.decoy_tag.clone();
    let generate_decoys = params.generate_decoys;
    let peptides = params.digest(&fasta);

    log::info!("digested {} target/decoy peptides", peptides.len());

    let bytes = sage_cloudpath::parquet::serialize_peptides(&peptides, &decoy_tag, generate_decoys)?;

    // `sage_cloudpath::to_url` canonicalizes the whole path, which requires
    // it to already exist — fine for inputs, wrong for a file we're about
    // to create. Canonicalize the (existing) parent directory instead, same
    // pattern as `Runner::make_path`'s `output_directory.join(...)`.
    let output_path = std::path::Path::new(output);
    let parent = match output_path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => std::path::Path::new("."),
    };
    std::fs::create_dir_all(parent)?;
    let file_name = output_path
        .file_name()
        .and_then(|n| n.to_str())
        .ok_or_else(|| anyhow::anyhow!("invalid output path: {output}"))?;
    let url = sage_cloudpath::Url::from_directory_path(parent.canonicalize()?)
        .map_err(|_| anyhow::anyhow!("invalid output directory: {}", parent.display()))?
        .join(file_name)?;
    sage_cloudpath::write_bytes_sync(&url, bytes)?;

    Ok(())
}
