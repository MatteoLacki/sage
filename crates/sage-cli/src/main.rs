use clap::{value_parser, Arg, Command, ValueHint};
use rayon::ThreadPoolBuilder;
use sage_cli::input::Input;
use sage_cli::runner::Runner;

fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("sage")
        .version(clap::crate_version!())
        .author("Michael Lazear <michaellazear92@gmail.com>")
        .about("\u{1F52E} Sage \u{1F9D9} - Proteomics searching so fast it feels like magic!")
        .arg(
            Arg::new("parameters")
                .required(true)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help("Path to configuration parameters (JSON file)")
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("mzml_paths")
                .num_args(1..)
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Paths to mzML, MGF, or .pmsms directories to process. Overrides files listed in the \
                     configuration file. Not used together with --pmsms/--precursors.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("pmsms")
                .long("pmsms")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to a pmsms.mmappet dataset (with an `mz` column already materialized). \
                     Must be given together with --precursors, as an alternative to the \
                     positional spectra paths.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("precursors")
                .long("precursors")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to a precursors table, .parquet or .mmappet (paired with --pmsms).",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("predicted-rt")
                .long("predicted-rt")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to a parquet file of externally-predicted peptide RT \
                     (columns: sequence, rt). Requires `rt_tol_sec` to also be \
                     configured. Independent of --predicted-iim.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("predicted-iim")
                .long("predicted-iim")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to a parquet file of externally-predicted peptide IIM \
                     (columns: sequence, charge, iim). Requires `mobility_tol` to \
                     also be configured. Independent of --predicted-rt.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("predicted-fragment-intensity-index")
                .long("predicted-fragment-intensity-index")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to a job-scoped parquet index (columns: sequence, charge, \
                     start, end) into a fragment-intensity PredictionCache's shared \
                     arrays.mmappet (git/featureprediction's \
                     export_fragment_intensity_for_sage). Feature-only -- computes \
                     Feature::ms2_entropy_similarity, no hard filtering. Must be given \
                     together with --predicted-fragment-intensity-cache.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("predicted-fragment-intensity-cache")
                .long("predicted-fragment-intensity-cache")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to the shared arrays.mmappet directory that \
                     --predicted-fragment-intensity-index's (start, end) ranges \
                     address (a subdirectory of the fragment-intensity \
                     PredictionCache -- not the whole cache directory, and not its \
                     index.sqlite3).",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("unimod-db-path")
                .long("unimod-db-path")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to a Unimod modifications CSV (id,name,mono_mass), overriding the \
                     table embedded in the binary. Lets `static_mods`/`variable_mods` reference \
                     a modification as `\"UNIMOD:<id>\"` instead of a raw mass delta.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("fasta")
                .short('f')
                .long("fasta")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path to FASTA database. Overrides the FASTA file \
                     specified in the configuration file.",
                )
                .value_hint(ValueHint::FilePath),
        )
        .arg(
            Arg::new("output_directory")
                .short('o')
                .long("output_directory")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .help(
                    "Path where search and quant results will be written. \
                     Overrides the directory specified in the configuration file.",
                )
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("batch-size")
                .long("batch-size")
                .value_parser(value_parser!(u16).range(1..))
                .help("Number of files to load and search in parallel (default = # of CPUs/2)")
                .value_hint(ValueHint::Other),
        )
        .arg(
            Arg::new("parquet")
                .long("parquet")
                .action(clap::ArgAction::SetTrue)
                .help("Write search output in parquet format instead of tsv"),
        )
        .arg(
            Arg::new("annotate-matches")
                .long("annotate-matches")
                .action(clap::ArgAction::SetTrue)
                .help("Write matched fragments output file."),
        )
        .arg(
            Arg::new("write-pin")
                .long("write-pin")
                .action(clap::ArgAction::SetTrue)
                .help("Write percolator-compatible `.pin` output files"),
        )
        .arg(
            Arg::new("write-report")
                .long("write-report")
                .action(clap::ArgAction::SetTrue)
                .help("Write `.html` report file"),
        )
        .arg(
            Arg::new("disable-telemetry")
                .long("disable-telemetry-i-dont-want-to-improve-sage")
                .action(clap::ArgAction::SetFalse)
                .help("Disable sending telemetry data"),
        )
        .arg(
            Arg::new("stack-size")
                .long("stack-size")
                .value_parser(value_parser!(u32).range(1..))
                .help("Set Rayon worker thread stack size in MiB (default: 2 MiB)")
                .value_hint(ValueHint::Other),
        )
        .help_template(
            "{usage-heading} {usage}\n\n\
             {about-with-newline}\n\
             Written by {author-with-newline}Version {version}\n\n\
             {all-args}{after-help}",
        )
        .get_matches();

    let stack_size_mib = matches.get_one::<u32>("stack-size").copied().unwrap_or(2);
    let stack_size_bytes = stack_size_mib as usize * 1024 * 1024;
    log::trace!(
        "setting Rayon worker thread stack size to {} MiB",
        stack_size_mib
    );
    ThreadPoolBuilder::new()
        .stack_size(stack_size_bytes)
        .build_global()
        .expect("configure Rayon pool");

    let parallel = matches
        .get_one::<u16>("batch-size")
        .copied()
        .unwrap_or_else(|| num_cpus::get() as u16 / 2) as usize;

    let parquet = matches.get_one::<bool>("parquet").copied().unwrap_or(false);
    let send_telemetry = matches
        .get_one::<bool>("disable-telemetry")
        .copied()
        .unwrap_or(true);

    let input = Input::from_arguments(matches)?;

    let runner = input
        .build()
        .and_then(|parameters| Runner::new(parameters, parallel))?;

    let tel = runner.run(parallel, parquet)?;

    if send_telemetry {
        tel.send();
    }

    Ok(())
}
