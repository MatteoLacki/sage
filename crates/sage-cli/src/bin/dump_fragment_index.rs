use clap::{Arg, ArgAction, Command, ValueHint};
use sage_core::database::Builder;

/// Build the real `IndexedDatabase` for a FASTA database and dump a 2D
/// occupancy histogram of its (peptide monoisotopic mass, fragment neutral
/// mass) pairs, binned geometrically so every bin spans a fixed ppm width.
///
/// Takes the same digestion-only `database` config section `dump_peptides`
/// does (enzyme, static/variable mods, mass bounds, decoy_tag) plus the
/// keys the fragment set itself depends on (`ion_kinds`, `min_ion_index`).
/// `bucket_size` only orders fragments within a page, which the histogram
/// cannot see.
///
/// The histogram is accumulated here rather than dumped as ~150M raw pairs:
/// a true 5ppm x 5ppm grid over the full mass ranges is ~3.2e11 cells, so
/// bins are folded by an integer factor into a grid of about
/// `--target-pixels` on a side. The integer fold keeps every pixel edge on
/// an exact ppm-bin boundary.
///
/// Output is an mmappet directory in the single-flat-column convention
/// `timstofu`'s `binary/array_serialization.py` reads: `schema.txt`
/// (`uint32 x`), `0.bin` (C-order counts, precursor-major), `shape.txt`
/// (`<nx>x<ny>`), plus a `grid.json` sidecar holding the bin parameters
/// needed to reconstruct mass axes.
fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Info)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("dump_fragment_index")
        .version(clap::crate_version!())
        .about(
            "Dump a ppm-binned precursor-mass x fragment-mass histogram of sage's fragment index",
        )
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
            Arg::new("output")
                .short('o')
                .long("output")
                .value_parser(clap::builder::NonEmptyStringValueParser::new())
                .default_value("fragment_index_grid.mmappet")
                .help("Output mmappet directory path")
                .value_hint(ValueHint::DirPath),
        )
        .arg(
            Arg::new("ppm")
                .long("ppm")
                .value_parser(clap::value_parser!(f64))
                .default_value("5.0")
                .help("Width of one geometric mass bin, in ppm, on both axes"),
        )
        .arg(
            Arg::new("bin-da")
                .long("bin-da")
                .value_parser(clap::value_parser!(f64))
                .help(
                    "Bin linearly with this width in Da on both axes, instead of geometrically \
                     by --ppm",
                ),
        )
        .arg(
            Arg::new("target-pixels")
                .long("target-pixels")
                .value_parser(clap::value_parser!(u64).range(1..))
                .default_value("8192")
                .help(
                    "Upper bound on the grid side length; ppm bins are folded by the smallest \
                     integer factor that fits each axis within it",
                ),
        )
        .arg(
            Arg::new("include-decoys")
                .long("include-decoys")
                .action(ArgAction::SetTrue)
                .help("Count decoy peptides too (default: target peptides only)"),
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
    let output = std::path::PathBuf::from(matches.get_one::<String>("output").unwrap());
    let ppm = *matches.get_one::<f64>("ppm").unwrap();
    let bin_da = matches.get_one::<f64>("bin-da").copied();
    let target_pixels = *matches.get_one::<u64>("target-pixels").unwrap();
    let include_decoys = matches.get_flag("include-decoys");

    if !(ppm > 0.0) {
        anyhow::bail!("--ppm must be positive, got {ppm}");
    }
    if let Some(width) = bin_da {
        if !(width > 0.0) {
            anyhow::bail!("--bin-da must be positive, got {width}");
        }
    }

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

    let db = params.build(fasta);
    log::info!(
        "indexed {} peptides into {} fragments",
        db.peptides.len(),
        db.size()
    );

    let keep_peptide: Vec<bool> = db
        .peptides
        .iter()
        .map(|peptide| include_decoys || !peptide.decoy)
        .collect();
    let precursor_mass: Vec<f64> = db
        .peptides
        .iter()
        .map(|peptide| peptide.monoisotopic as f64)
        .collect();

    let precursor_range = mass_range(
        precursor_mass
            .iter()
            .zip(&keep_peptide)
            .filter(|(_, keep)| **keep)
            .map(|(mass, _)| *mass),
    );
    let fragment_range = mass_range(
        db.fragment_mzs
            .iter()
            .zip(&db.fragment_peptide_indices)
            .filter(|(_, peptide_index)| keep_peptide[peptide_index.0 as usize])
            .map(|(fragment_mz, _)| *fragment_mz as f64),
    );

    let (precursor_min, precursor_max) =
        precursor_range.ok_or_else(|| anyhow::anyhow!("no peptides left to histogram"))?;
    let (fragment_min, fragment_max) =
        fragment_range.ok_or_else(|| anyhow::anyhow!("no fragments left to histogram"))?;

    let binning = match bin_da {
        Some(width) => Binning::Linear { width },
        None => Binning::Geometric {
            ppm,
            log_step: (1.0 + ppm * 1e-6).ln(),
        },
    };
    let precursor_axis = Axis::new(precursor_min, precursor_max, binning, target_pixels);
    let fragment_axis = Axis::new(fragment_min, fragment_max, binning, target_pixels);
    log::info!(
        "precursor axis: {:.4}-{:.4} Da, m0 {}, {} bins, fold {}, {} pixels",
        precursor_min,
        precursor_max,
        precursor_axis.m0,
        precursor_axis.n_bins,
        precursor_axis.fold,
        precursor_axis.n_pixels
    );
    log::info!(
        "fragment axis: {:.4}-{:.4} Da, m0 {}, {} bins, fold {}, {} pixels",
        fragment_min,
        fragment_max,
        fragment_axis.m0,
        fragment_axis.n_bins,
        fragment_axis.fold,
        fragment_axis.n_pixels
    );

    let precursor_pixel: Vec<u32> = precursor_mass
        .iter()
        .map(|mass| precursor_axis.pixel(*mass) as u32)
        .collect();

    let mut counts = vec![0u32; (precursor_axis.n_pixels * fragment_axis.n_pixels) as usize];
    let mut counted = 0u64;
    for (fragment_mz, peptide_index) in db.fragment_mzs.iter().zip(&db.fragment_peptide_indices) {
        let peptide_index = peptide_index.0 as usize;
        if !keep_peptide[peptide_index] {
            continue;
        }
        let row = precursor_pixel[peptide_index] as u64;
        let column = fragment_axis.pixel(*fragment_mz as f64);
        counts[(row * fragment_axis.n_pixels + column) as usize] += 1;
        counted += 1;
    }
    let occupied = counts.iter().filter(|count| **count > 0).count();
    log::info!(
        "histogrammed {counted} (precursor, fragment) pairs into {occupied} occupied pixels"
    );

    write_grid(
        &output,
        &counts,
        &precursor_axis,
        &fragment_axis,
        include_decoys,
        counted,
    )?;
    log::info!("wrote {}", output.display());

    Ok(())
}

/// One mass axis: `fold` adjacent bins per output pixel, bins starting at
/// `m0`, the integer floor of the smallest mass, so bin edges are
/// reproducible from `grid.json` alone.
#[derive(Copy, Clone)]
enum Binning {
    /// Bin `k` spans `[m0 * r^k, m0 * r^(k+1))` for `r = 1 + ppm * 1e-6`:
    /// constant ppm width, pixel index linear in `log(mass)`.
    Geometric { ppm: f64, log_step: f64 },
    /// Bin `k` spans `[m0 + k*width, m0 + (k+1)*width)`: constant Da width,
    /// pixel index linear in mass.
    Linear { width: f64 },
}

#[derive(Copy, Clone)]
struct Axis {
    m0: f64,
    binning: Binning,
    n_bins: u64,
    fold: u64,
    n_pixels: u64,
}

impl Axis {
    fn new(min_mass: f64, max_mass: f64, binning: Binning, target_pixels: u64) -> Self {
        let m0 = min_mass.floor().max(1.0);
        let n_bins = match binning {
            Binning::Geometric { log_step, .. } => ((max_mass / m0).ln() / log_step).floor() as u64,
            Binning::Linear { width } => ((max_mass - m0) / width).floor() as u64,
        } + 1;
        let fold = n_bins.div_ceil(target_pixels).max(1);
        Self {
            m0,
            binning,
            n_bins,
            fold,
            n_pixels: n_bins.div_ceil(fold),
        }
    }

    fn pixel(&self, mass: f64) -> u64 {
        let bin = match self.binning {
            Binning::Geometric { log_step, .. } => (mass / self.m0).ln() / log_step,
            Binning::Linear { width } => (mass - self.m0) / width,
        };
        ((bin.floor().max(0.0) as u64) / self.fold).min(self.n_pixels - 1)
    }

    fn as_meta(&self) -> serde_json::Value {
        let mut meta = serde_json::json!({
            "m0": self.m0,
            "n_bins": self.n_bins,
            "fold": self.fold,
            "n_pixels": self.n_pixels,
        });
        match self.binning {
            Binning::Geometric { ppm, .. } => {
                meta["scale"] = "geometric".into();
                meta["ppm"] = ppm.into();
            }
            Binning::Linear { width } => {
                meta["scale"] = "linear".into();
                meta["bin_da"] = width.into();
            }
        }
        meta
    }
}

fn mass_range(masses: impl Iterator<Item = f64>) -> Option<(f64, f64)> {
    masses.fold(None, |range, mass| match range {
        None => Some((mass, mass)),
        Some((lo, hi)) => Some((lo.min(mass), hi.max(mass))),
    })
}

fn write_grid(
    output: &std::path::Path,
    counts: &[u32],
    precursor_axis: &Axis,
    fragment_axis: &Axis,
    include_decoys: bool,
    counted: u64,
) -> anyhow::Result<()> {
    use std::io::Write;

    std::fs::create_dir_all(output)?;
    // mmappet's own `schema_to_str` joins lines without a trailing newline;
    // match that byte-for-byte so `shape.txt`-style datasets written here
    // are indistinguishable from ones written by the Python writer.
    std::fs::write(output.join("schema.txt"), "uint32 x")?;
    std::fs::write(
        output.join("shape.txt"),
        format!("{}x{}", precursor_axis.n_pixels, fragment_axis.n_pixels),
    )?;

    let mut column = std::io::BufWriter::new(std::fs::File::create(output.join("0.bin"))?);
    for count in counts {
        column.write_all(&count.to_le_bytes())?;
    }
    column.flush()?;

    let grid = serde_json::json!({
        "title": "SAGE fragment index occupancy",
        "x_label": "peptide monoisotopic mass [Da]",
        "y_label": "fragment neutral mass [Da]",
        "targets_only": !include_decoys,
        "n_pairs": counted,
        "precursor": precursor_axis.as_meta(),
        "fragment": fragment_axis.as_meta(),
    });
    std::fs::write(
        output.join("grid.json"),
        serde_json::to_string_pretty(&grid)?,
    )?;
    Ok(())
}
