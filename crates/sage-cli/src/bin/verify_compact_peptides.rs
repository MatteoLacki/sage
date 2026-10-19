use clap::{Arg, ArgAction, Command, ValueHint};
use rayon::prelude::*;
use sage_cli::process_memory::{peak_rss_bytes, reset_peak_rss};
use sage_core::compact::{drop_in_parallel, CompactPeptide, CompactPeptides};
use sage_core::database::Builder;
use sage_core::enzyme::Position;
use sage_core::peptide::Peptide;
use std::time::Instant;

/// Check `CompactPeptides` against sage's own digest, at full scale.
///
/// Per peptide length: digest with sage, compress, and require every row to
/// materialize back into the identical `Peptide` (all fields but `proteins`), and
/// every `--protein-check-every`-th row to resolve the identical protein list.
/// Then merge the lengths and require sage's row order between adjacent rows.
///
/// With `--full`, also digest the whole database at once (it must fit in memory)
/// and require the merged table to equal it row for row, proteins included. A
/// peptide merged from several digests keeps an arbitrary duplicate's `position`,
/// `semi_enzymatic` and `missed_cleavages`, so those are counted, not failed.
///
/// With `--direct` (non-specific digests), first build the table with
/// `CompactPeptides::build_nonspecific`, then require it to equal the per-length
/// table row for row: mass bits, length, decoy flag, variable sites and the residues
/// at `start`. Differing `start` offsets and `position`s are counted, not failed.
///
/// With `--views N`, only build the table (`CompactPeptides::build`) and require
/// every N-th row's `CompactView` to read as its materialised `Peptide` does
/// (`peptide_view::first_difference`: residues, modification bits, termini, mass,
/// decoy, missed cleavages, every ion).
///
/// Prints a TSV per length to stdout and `#`-prefixed summary lines; exits
/// non-zero on any mismatch.
fn main() -> anyhow::Result<()> {
    env_logger::Builder::default()
        .filter_level(log::LevelFilter::Error)
        .parse_env(env_logger::Env::default().filter_or("SAGE_LOG", "error,sage=info"))
        .init();

    let matches = Command::new("verify_compact_peptides")
        .version(clap::crate_version!())
        .about("Check compact peptide rows against sage's own digest")
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
        .arg(
            Arg::new("protein-check-every")
                .long("protein-check-every")
                .value_parser(clap::value_parser!(usize))
                .default_value("1000")
                .help("Resolve and compare the proteins of every N-th row of each length"),
        )
        .arg(
            Arg::new("full")
                .long("full")
                .action(ArgAction::SetTrue)
                .help("Also compare against the whole database digested at once"),
        )
        .arg(
            Arg::new("direct")
                .long("direct")
                .action(ArgAction::SetTrue)
                .help("Also build the table directly (non-specific digests) and compare"),
        )
        .arg(
            Arg::new("views")
                .long("views")
                .value_parser(clap::value_parser!(usize))
                .help("Only compare every N-th row's view against its materialised peptide"),
        )
        .get_matches();

    let fasta_path = matches.get_one::<String>("fasta").unwrap();
    let protein_check_every = (*matches.get_one::<usize>("protein-check-every").unwrap()).max(1);
    let full = matches.get_flag("full");
    let direct = matches.get_flag("direct");

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

    if let Some(&every) = matches.get_one::<usize>("views") {
        return check_views(&params, &fasta, every.max(1));
    }

    let direct_table = match direct {
        true => {
            reset_peak_rss()?;
            let start = Instant::now();
            let direct_table = CompactPeptides::build_nonspecific(&params, &fasta)?;
            println!(
                "# direct build: {} rows ({:.3} GB) in {:.2} s, {} adjacent rows out of sage's order, peak rss {:.3} GB",
                direct_table.rows.len(),
                direct_table.rows.bytes() as f64 / 1e9,
                start.elapsed().as_secs_f64(),
                direct_table.row_order_violations(),
                peak_rss_bytes()? as f64 / 1e9,
            );
            Some(direct_table)
        }
        false => None,
    };

    let mut table = CompactPeptides::new(&params, &fasta)?;
    println!(
        "length\trows\tcompact_gb\tfield_mismatches\tproteins_checked\tprotein_mismatches\tdigest_s\tcompress_s\tcheck_s\tpeak_rss_gb"
    );
    let mut failures = 0;
    let mut runs = Vec::new();
    for length in table.lengths() {
        reset_peak_rss()?;
        let digest_start = Instant::now();
        let peptides = table.digest_length(&fasta, length);
        let digest_seconds = digest_start.elapsed().as_secs_f64();

        let compress_start = Instant::now();
        let rows = table.compress(&peptides)?;
        let compress_seconds = compress_start.elapsed().as_secs_f64();

        let check_start = Instant::now();
        let field_mismatches = rows
            .par_iter()
            .zip(&peptides)
            .filter(|(row, peptide)| !same_fields(&table.materialize(row), peptide))
            .count();
        let checked_rows: Vec<usize> = (0..rows.len()).step_by(protein_check_every).collect();
        let protein_mismatches = protein_mismatches(&table, &rows, &peptides, &checked_rows);
        let check_seconds = check_start.elapsed().as_secs_f64();

        println!(
            "{}\t{}\t{:.3}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.3}",
            length,
            rows.len(),
            compact_gb(rows.len()),
            field_mismatches,
            checked_rows.len(),
            protein_mismatches,
            digest_seconds,
            compress_seconds,
            check_seconds,
            peak_rss_bytes()? as f64 / 1e9,
        );
        failures += field_mismatches + protein_mismatches;
        runs.push(rows);
        drop_in_parallel(peptides);
    }

    reset_peak_rss()?;
    let merge_start = Instant::now();
    table.merge_length_runs(runs)?;
    let order_violations = table.row_order_violations();
    println!(
        "# merged {} rows ({:.3} GB) in {:.2} s, {} adjacent rows out of sage's order, peak rss {:.3} GB",
        table.rows.len(),
        table.rows.bytes() as f64 / 1e9,
        merge_start.elapsed().as_secs_f64(),
        order_violations,
        peak_rss_bytes()? as f64 / 1e9,
    );
    failures += order_violations;

    if let Some(direct_table) = &direct_table {
        failures += direct_table.row_order_violations();
        failures += compare_with_direct_build(&table, direct_table);
    }

    if full {
        failures += compare_with_full_digest(&table, &params.digest(&fasta));
    }

    println!("# {}", if failures == 0 { "PASS" } else { "FAIL" });
    if failures > 0 {
        std::process::exit(1);
    }
    Ok(())
}

fn compact_gb(rows: usize) -> f64 {
    (rows * std::mem::size_of::<CompactPeptide>()) as f64 / 1e9
}

/// Every field but `proteins`, which `materialize` leaves to `resolve_proteins`.
fn same_fields(materialized: &Peptide, expected: &Peptide) -> bool {
    Peptide {
        proteins: Vec::new(),
        ..expected.clone()
    } == *materialized
        && materialized.modifications == expected.modifications
}

fn protein_mismatches(
    table: &CompactPeptides,
    rows: &[CompactPeptide],
    expected: &[Peptide],
    checked_rows: &[usize],
) -> usize {
    let mut resolved: Vec<Peptide> = checked_rows
        .par_iter()
        .map(|&row| table.materialize(&rows[row]))
        .collect();
    table.resolve_proteins(&mut resolved);
    resolved
        .iter()
        .zip(checked_rows)
        .filter(|(peptide, &row)| peptide.proteins != expected[row].proteins)
        .count()
}

/// Row-for-row comparison with the whole database digested at once. Returns the
/// failure count; arbitrary duplicate attributes are reported but not failures.
fn compare_with_full_digest(table: &CompactPeptides, expected: &[Peptide]) -> usize {
    if table.rows.len() != expected.len() {
        println!(
            "# full digest: {} rows, compact {} rows",
            expected.len(),
            table.rows.len()
        );
        return 1;
    }
    let start = Instant::now();
    let mut materialized: Vec<Peptide> = (0..table.rows.len())
        .into_par_iter()
        .map(|row| table.materialize(&table.rows.get(row)))
        .collect();
    table.resolve_proteins(&mut materialized);

    let field_mismatches = materialized
        .par_iter()
        .zip(expected)
        .filter(|(actual, expected)| {
            without_duplicate_attributes(actual) != without_duplicate_attributes(expected)
                || actual.modifications != expected.modifications
        })
        .count();
    let (position, semi_enzymatic, missed_cleavages) = materialized
        .par_iter()
        .zip(expected)
        .map(|(actual, expected)| {
            (
                (actual.position != expected.position) as usize,
                (actual.semi_enzymatic != expected.semi_enzymatic) as usize,
                (actual.missed_cleavages != expected.missed_cleavages) as usize,
            )
        })
        .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
    println!(
        "# full digest: {} rows, {} mismatches in everything but duplicate attributes \
         (incl. proteins); duplicate attributes differing: position {}, semi_enzymatic {}, \
         missed_cleavages {}; {:.2} s",
        expected.len(),
        field_mismatches,
        position,
        semi_enzymatic,
        missed_cleavages,
        start.elapsed().as_secs_f64(),
    );
    field_mismatches
}

/// Row-for-row comparison of the per-length table with the direct build. Returns
/// the failure count; `start` (same residues, another occurrence) and `position`
/// differences are reported, not failures.
fn compare_with_direct_build(by_length: &CompactPeptides, direct: &CompactPeptides) -> usize {
    if by_length.rows.len() != direct.rows.len() {
        println!(
            "# direct build: {} rows, per-length {} rows",
            direct.rows.len(),
            by_length.rows.len()
        );
        return 1;
    }
    let (mismatches, other_starts, other_positions) = (0..direct.rows.len())
        .into_par_iter()
        .map(|row| {
            let (a, b) = (by_length.rows.get(row), direct.rows.get(row));
            let same = a.monoisotopic.to_bits() == b.monoisotopic.to_bits()
                && by_length.same_row_apart_from_occurrence(&a, direct, &b);
            (
                !same as usize,
                (a.start() != b.start()) as usize,
                (a.position() != b.position()) as usize,
            )
        })
        .reduce(|| (0, 0, 0), |a, b| (a.0 + b.0, a.1 + b.1, a.2 + b.2));
    println!(
        "# direct build vs per-length: {} rows, {} mismatches; start at another occurrence {}, position differing {}",
        direct.rows.len(),
        mismatches,
        other_starts,
        other_positions,
    );
    mismatches
}

/// The attributes a peptide merged from several digests takes from whichever
/// duplicate survives `reorder_peptides`' unstable sort.
fn without_duplicate_attributes(peptide: &Peptide) -> Peptide {
    Peptide {
        position: Position::Internal,
        semi_enzymatic: false,
        missed_cleavages: 0,
        ..peptide.clone()
    }
}

/// Every `every`-th row of the table: its view against its materialised peptide.
fn check_views(
    params: &sage_core::database::Parameters,
    fasta: &sage_core::fasta::Fasta,
    every: usize,
) -> anyhow::Result<()> {
    let start = Instant::now();
    let table = CompactPeptides::build(params, fasta)?;
    let built = start.elapsed().as_secs_f64();
    let start = Instant::now();
    let rows: Vec<usize> = (0..table.rows.len()).step_by(every).collect();
    let differences: Vec<(usize, String)> = rows
        .par_iter()
        .filter_map(|&row| {
            let peptide = table.materialize(&table.rows.get(row));
            sage_core::peptide_view::first_difference(&table.view(row), &peptide)
                .map(|difference| (row, format!("{peptide}: {difference}")))
        })
        .collect();
    for (row, difference) in differences.iter().take(20) {
        println!("# row {row}: {difference}");
    }
    println!(
        "# views: {} of {} rows checked (every {every}), {} differ; table built in {built:.1} s, checked in {:.1} s",
        rows.len(),
        table.rows.len(),
        differences.len(),
        start.elapsed().as_secs_f64()
    );
    anyhow::ensure!(differences.is_empty(), "{} views differ", differences.len());
    Ok(())
}
