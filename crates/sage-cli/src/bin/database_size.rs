use clap::{Arg, Command, ValueHint};
use sage_core::database::{Builder, IndexedDatabase, Parameters};
use sage_core::enzyme::EnzymeParameters;
use sage_core::peptide::Peptide;
use std::sync::Arc;
use std::time::Instant;

/// Measure how big sage's `IndexedDatabase` gets for a FASTA plus a `database`
/// config section, by running sage's own `digest` and `build_from_peptides`
/// one peptide length at a time.
///
/// Per-length builds sum to the exact full-database counts: dedup only ever
/// merges peptides with identical sequences, and a decoy is only dropped when
/// it equals a target sequence, so peptides of different lengths never
/// interact. Peak memory is therefore that of the largest single length, not
/// of the whole database.
///
/// Prints a TSV to stdout, one row per length plus a `total` row. Byte columns
/// are estimates from struct sizes and allocation capacities, ignoring
/// allocator rounding. `peak_rss_gb` is the process high-water mark during
/// that length's digest+build (reset between lengths via
/// `/proc/self/clear_refs`), so it includes the FASTA held throughout.
fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("database_size")
        .version(clap::crate_version!())
        .about("Report peptide/fragment counts and memory of sage's index, per peptide length")
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
                    "Optional path to a JSON `database` section (same shape as `sage`'s own \
                     config `database` key) to override defaults. Omit to use sage's defaults.",
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
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();

    // Must happen before parsing --config below: overrides the embedded
    // default table `static_mods`/`variable_mods`'s `UNIMOD:<id>`
    // references get resolved against.
    if let Some(unimod_db_path) = matches.get_one::<String>("unimod-db-path") {
        sage_core::unimod::set_active_table_from_path(std::path::Path::new(unimod_db_path))
            .map_err(|e| anyhow::anyhow!(e))?;
    }

    let mut builder: Builder = match matches.get_one::<String>("config") {
        Some(path) => {
            let mut value: serde_json::Value = serde_json::from_slice(&std::fs::read(path)?)?;
            sage_cli::input::resolve_unimod_refs_in_database(&mut value)?;
            serde_json::from_value(value)?
        }
        None => Builder::default(),
    };
    builder.fasta = Some(fasta_path.clone());
    let params = builder.make_parameters();

    let fasta =
        sage_cloudpath::util::read_fasta(&params.fasta, &params.decoy_tag, params.generate_decoys)?;

    let enzyme: EnzymeParameters = params.enzyme.clone().into();

    println!(
        "length\ttarget_peptides\tdecoy_peptides\tfragments\tpages\tfragment_index_gb\tpeptides_gb\tdigest_s\tbuild_s\tpeak_rss_gb"
    );
    let mut total = LengthSize::default();
    for length in enzyme.min_len..=enzyme.max_len {
        reset_peak_rss()?;
        let length_params = params_for_single_length(&params, length);

        let digest_start = Instant::now();
        let peptides = length_params.digest(&fasta);
        let digest_seconds = digest_start.elapsed().as_secs_f64();

        let build_start = Instant::now();
        let db = length_params.build_from_peptides(peptides);
        let build_seconds = build_start.elapsed().as_secs_f64();

        let size = LengthSize::measure(&db, digest_seconds, build_seconds, peak_rss_bytes()?);
        size.print(&length.to_string());
        total.accumulate(&size);
    }
    total.print("total");

    Ok(())
}

fn params_for_single_length(params: &Parameters, length: usize) -> Parameters {
    let mut length_params = params.clone();
    length_params.enzyme.min_len = Some(length);
    length_params.enzyme.max_len = Some(length);
    length_params
}

#[derive(Default)]
struct LengthSize {
    target_peptides: usize,
    decoy_peptides: usize,
    fragments: usize,
    pages: usize,
    fragment_index_bytes: usize,
    peptides_bytes: usize,
    digest_seconds: f64,
    build_seconds: f64,
    peak_rss_bytes: usize,
}

impl LengthSize {
    fn measure(
        db: &IndexedDatabase,
        digest_seconds: f64,
        build_seconds: f64,
        peak_rss_bytes: usize,
    ) -> Self {
        let decoy_peptides = db.peptides.iter().filter(|peptide| peptide.decoy).count();
        let fragment_index_bytes = db.fragment_mzs.capacity() * std::mem::size_of::<f32>()
            + db.fragment_peptide_indices.capacity() * std::mem::size_of::<u32>()
            + db.min_value.capacity() * std::mem::size_of::<f32>();
        let peptides_bytes = db.peptides.capacity() * std::mem::size_of::<Peptide>()
            + db.peptides.iter().map(peptide_heap_bytes).sum::<usize>();
        Self {
            target_peptides: db.peptides.len() - decoy_peptides,
            decoy_peptides,
            fragments: db.size(),
            pages: db.min_value.len(),
            fragment_index_bytes,
            peptides_bytes,
            digest_seconds,
            build_seconds,
            peak_rss_bytes,
        }
    }

    fn accumulate(&mut self, other: &Self) {
        self.target_peptides += other.target_peptides;
        self.decoy_peptides += other.decoy_peptides;
        self.fragments += other.fragments;
        self.pages += other.pages;
        self.fragment_index_bytes += other.fragment_index_bytes;
        self.peptides_bytes += other.peptides_bytes;
        self.digest_seconds += other.digest_seconds;
        self.build_seconds += other.build_seconds;
        self.peak_rss_bytes = self.peak_rss_bytes.max(other.peak_rss_bytes);
    }

    fn print(&self, label: &str) {
        const GB: f64 = 1e9;
        println!(
            "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.2}\t{:.2}\t{:.3}",
            label,
            self.target_peptides,
            self.decoy_peptides,
            self.fragments,
            self.pages,
            self.fragment_index_bytes as f64 / GB,
            self.peptides_bytes as f64 / GB,
            self.digest_seconds,
            self.build_seconds,
            self.peak_rss_bytes as f64 / GB,
        );
    }
}

/// Heap owned by one `Peptide`: its `Arc<[u8]>` sequence (two refcounts plus
/// bytes), its modification masses, and its protein-name pointers. The
/// protein names themselves are shared with the FASTA and not counted.
fn peptide_heap_bytes(peptide: &Peptide) -> usize {
    const ARC_REFCOUNTS_BYTES: usize = 2 * std::mem::size_of::<usize>();
    ARC_REFCOUNTS_BYTES
        + peptide.sequence.len()
        + peptide.modifications.capacity() * std::mem::size_of::<f32>()
        + peptide.proteins.capacity() * std::mem::size_of::<Arc<str>>()
}

/// Linux: writing "5" to `/proc/self/clear_refs` resets `VmHWM` to the current RSS.
fn reset_peak_rss() -> std::io::Result<()> {
    std::fs::write("/proc/self/clear_refs", "5")
}

fn peak_rss_bytes() -> anyhow::Result<usize> {
    let status = std::fs::read_to_string("/proc/self/status")?;
    let line = status
        .lines()
        .find(|line| line.starts_with("VmHWM:"))
        .ok_or_else(|| anyhow::anyhow!("no VmHWM line in /proc/self/status"))?;
    let kilobytes: usize = line
        .trim_start_matches("VmHWM:")
        .trim()
        .trim_end_matches("kB")
        .trim()
        .parse()?;
    Ok(kilobytes * 1024)
}
