use anyhow::{ensure, Context};
use clap::ArgMatches;
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::util::PmsmsPaths;
use sage_cloudpath::Url;
use sage_core::scoring::ScoreType;
use sage_core::spline::{FragmentTolSpline, ValueTolSpline};
use sage_core::{
    database::{Builder, Parameters},
    lfq::LfqSettings,
    mass::Tolerance,
    tmt::Isobaric,
};
use serde::{Deserialize, Serialize};

#[derive(Serialize, Clone)]
/// Actual search parameters - may include overrides or default values not set by user
pub struct Search {
    pub version: String,
    pub database: Parameters,
    pub quant: QuantSettings,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub fragment_tol_spline: Option<FragmentTolSpline>,
    /// Predicted-RT/IIM candidate filtering — see `predicted_properties` doc
    /// on [`Input`]. `rt_tol`'s spline *values* here are already converted
    /// to minutes; its grid stays in `scan_start_time`'s native minutes too.
    pub predicted_properties: Option<String>,
    pub rt_tol: Option<ValueTolSpline>,
    pub mobility_tol: Option<ValueTolSpline>,
    pub precursor_charge: (u8, u8),
    pub override_precursor_charge: bool,
    pub isotope_errors: (i8, i8),
    pub deisotope: bool,
    pub chimera: bool,
    pub wide_window: bool,
    pub min_peaks: usize,
    pub max_peaks: usize,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: u16,
    pub report_psms: usize,
    pub predict_rt: bool,
    pub mzml_paths: Vec<Url>,
    pub pmsms_paths: Option<PmsmsPaths>,
    pub output_paths: Vec<Url>,
    pub bruker_config: BrukerProcessingConfig,
    pub protein_grouping: bool,
    pub protein_grouping_peptide_fdr: f32,

    #[serde(skip_serializing)]
    pub output_directory: Url,

    #[serde(skip_serializing)]
    pub write_pin: bool,

    #[serde(skip_serializing)]
    pub write_report: bool,

    #[serde(skip_serializing)]
    pub annotate_matches: bool,

    pub score_type: ScoreType,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    pub database: Builder,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
    pub fragment_tol_spline: Option<FragmentTolSpline>,
    pub report_psms: Option<usize>,
    pub chimera: Option<bool>,
    pub wide_window: Option<bool>,
    pub min_peaks: Option<usize>,
    pub max_peaks: Option<usize>,
    pub max_fragment_charge: Option<u8>,
    pub min_matched_peaks: Option<u16>,
    pub precursor_charge: Option<(u8, u8)>,
    pub override_precursor_charge: Option<bool>,
    pub isotope_errors: Option<(i8, i8)>,
    pub deisotope: Option<bool>,
    pub quant: Option<QuantOptions>,
    pub predict_rt: Option<bool>,
    pub output_directory: Option<String>,
    pub mzml_paths: Option<Vec<String>>,
    pub pmsms: Option<String>,
    pub precursors: Option<String>,
    /// Path to a parquet file of externally-predicted peptide RT/IIM
    /// (columns: sequence, charge, rt, iim), used to reject candidates
    /// whose predicted RT/IIM falls outside `rt_tol_sec`/`mobility_tol` of
    /// the observed spectrum's values. Requires both to be set.
    pub predicted_properties: Option<String>,
    /// RT tolerance window as a function of observed `scan_start_time`
    /// (minutes — the spline's own grid is in that same unit, no
    /// conversion). Spline *values* (the window width) are in **seconds**
    /// (converted to minutes internally at load time) — unit spelled out
    /// in the name since it differs from SAGE's internal minutes. A flat
    /// (value-independent) window is just a 2-node spline with identical
    /// values at both nodes — there is no separate flat-tolerance type.
    pub rt_tol_sec: Option<ValueTolSpline>,
    /// IIM tolerance window as a function of observed
    /// `Precursor::inverse_ion_mobility` (1/K0 — unambiguous, no unit
    /// suffix needed, same unit for both the spline's grid and its values).
    pub mobility_tol: Option<ValueTolSpline>,
    pub bruker_config: Option<BrukerProcessingConfig>,
    pub protein_grouping: Option<bool>,
    pub protein_grouping_peptide_fdr: Option<f32>,

    pub annotate_matches: Option<bool>,
    pub write_pin: Option<bool>,
    pub write_report: Option<bool>,
    pub score_type: Option<ScoreType>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct LfqOptions {
    pub peak_scoring: Option<sage_core::lfq::PeakScoringStrategy>,
    pub integration: Option<sage_core::lfq::IntegrationStrategy>,
    pub spectral_angle: Option<f64>,
    pub ppm_tolerance: Option<f32>,
    pub mobility_pct_tolerance: Option<f32>,
    pub combine_charge_states: Option<bool>,
    pub peptide_q_value: Option<f32>,
}

impl From<LfqOptions> for LfqSettings {
    fn from(value: LfqOptions) -> LfqSettings {
        let default = LfqSettings::default();
        let settings = LfqSettings {
            peak_scoring: value.peak_scoring.unwrap_or(default.peak_scoring),
            integration: value.integration.unwrap_or(default.integration),
            spectral_angle: value.spectral_angle.unwrap_or(default.spectral_angle).abs(),
            ppm_tolerance: value.ppm_tolerance.unwrap_or(default.ppm_tolerance).abs(),
            peptide_q_value: value.peptide_q_value.unwrap_or(default.peptide_q_value),
            mobility_pct_tolerance: value
                .mobility_pct_tolerance
                .unwrap_or(default.mobility_pct_tolerance),
            combine_charge_states: value
                .combine_charge_states
                .unwrap_or(default.combine_charge_states),
        };
        if settings.ppm_tolerance > 20.0 {
            log::warn!("lfq_settings.ppm_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance > 4.0 {
            log::warn!("lfq_settings.mobility_pct_tolerance is higher than expected");
        }
        if settings.mobility_pct_tolerance < 0.05 {
            log::warn!("lfq_settings.mobility_pct_tolerance is smaller than expected");
        }
        if settings.spectral_angle < 0.50 {
            log::warn!("lfq_settings.spectral_angle is lower than expected");
        }
        if settings.peptide_q_value > 0.01 {
            log::info!("lfq_settings.peptide_q_value is higher than expected, expect increased runtime and memory usage");
        }
        if settings.peptide_q_value < 0.01 {
            log::warn!("lfq_settings.peptide_q_value is lower than expected, not all identified peptides will have MS1 intensities extracted");
        }

        settings
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct TmtOptions {
    pub level: Option<u8>,
    pub sn: Option<bool>,
}

#[derive(Copy, Clone, Serialize, Debug)]
pub struct TmtSettings {
    pub level: u8,
    pub sn: bool,
}

impl From<TmtOptions> for TmtSettings {
    fn from(value: TmtOptions) -> Self {
        let default = Self::default();
        Self {
            level: value.level.unwrap_or(default.level),
            sn: value.sn.unwrap_or(default.sn),
        }
    }
}

impl Default for TmtSettings {
    fn default() -> Self {
        Self {
            level: 3,
            sn: false,
        }
    }
}

#[derive(Serialize, Deserialize, Default, Debug)]
pub struct QuantOptions {
    pub tmt: Option<Isobaric>,
    #[serde(rename = "tmt_settings")]
    pub tmt_options: Option<TmtOptions>,

    pub lfq: Option<bool>,
    #[serde(rename = "lfq_settings")]
    pub lfq_options: Option<LfqOptions>,
}

#[derive(Serialize, Default, Clone)]
pub struct QuantSettings {
    pub tmt: Option<Isobaric>,
    pub tmt_settings: TmtSettings,
    pub lfq: bool,
    pub lfq_settings: LfqSettings,
}

impl From<QuantOptions> for QuantSettings {
    fn from(value: QuantOptions) -> Self {
        Self {
            tmt: value.tmt,
            tmt_settings: value.tmt_options.map(Into::into).unwrap_or_default(),

            lfq: value.lfq.unwrap_or(false),
            lfq_settings: value.lfq_options.map(Into::into).unwrap_or_default(),
        }
    }
}

/// Resolve one `database.static_mods`/`variable_mods` leaf value in place.
/// A `"UNIMOD:<id>"` string is replaced with the active table's exact
/// float, and `id` is recorded into `reverse` for `Peptide::Display`'s
/// output round-trip. A plain number is checked against the active table
/// too, but only to *reject* it if it coincides with a real entry --
/// silently letting a rounded/imprecise float through would mean a
/// modification that plainly *is* e.g. Carbamidomethyl never round-trips
/// as `[UNIMOD:4]` on output, an inconsistency the caller should fix by
/// writing `"UNIMOD:4"` explicitly instead.
fn resolve_mod_leaf(
    v: &mut serde_json::Value,
    reverse: &mut std::collections::HashMap<u32, u32>,
) -> anyhow::Result<()> {
    if let Some(s) = v.as_str() {
        let id_str = s.strip_prefix("UNIMOD:").ok_or_else(|| {
            anyhow::anyhow!(
                "modification `{s}` is a string but not `UNIMOD:<id>` -- \
                 mods must be a plain number or a `UNIMOD:<id>` reference"
            )
        })?;
        let id: u32 = id_str
            .parse()
            .map_err(|_| anyhow::anyhow!("invalid UNIMOD id in `{s}`"))?;
        let mass = sage_core::unimod::resolve(id).ok_or_else(|| {
            anyhow::anyhow!("UNIMOD:{id} not found in the active unimod table")
        })?;
        reverse.insert(mass.to_bits(), id);
        *v = serde_json::Value::from(mass);
        return Ok(());
    }

    if let Some(mass) = v.as_f64() {
        let mass = mass as f32;
        if let Some((id, canonical)) = sage_core::unimod::find_coincidental_match(mass) {
            anyhow::bail!(
                "modification value {mass} coincides with UNIMOD:{id} (mono_mass={canonical}) -- \
                 use \"UNIMOD:{id}\" instead of a raw float, so this modification round-trips \
                 correctly on output and matches what predictors that read UNIMOD notation expect"
            );
        }
    }

    Ok(())
}

/// Walks `static_mods`/`variable_mods` in a `database`-shaped raw config
/// JSON value (before it's deserialized into `Builder`), resolving
/// `"UNIMOD:<id>"` references and validating plain floats against the
/// active unimod table (embedded default, or `--unimod-db-path`'s
/// override -- callers must set that override, if any, before calling
/// this). `Builder.static_mods`/`variable_mods` keep their existing
/// `HashMap<String, f32>`/`HashMap<String, Vec<f32>>` shape unchanged --
/// this only ever produces a plain float for serde to deserialize into
/// them, same as if the config had been written that way directly.
///
/// `pub`, not `pub(crate)`: the `dump_peptides` binary (`sage-cli/src/
/// bin/dump_peptides.rs`) deserializes a `database`-shaped config
/// directly into `Builder` with no enclosing `Input`, and needs this same
/// resolution -- a separate `[[bin]]` target depends on `sage_cli` as an
/// external crate, so `pub(crate)` wouldn't be visible to it.
pub fn resolve_unimod_refs_in_database(database: &mut serde_json::Value) -> anyhow::Result<()> {
    let mut reverse = std::collections::HashMap::new();

    if let Some(static_mods) = database
        .get_mut("static_mods")
        .and_then(|v| v.as_object_mut())
    {
        for (_residue, v) in static_mods.iter_mut() {
            resolve_mod_leaf(v, &mut reverse)?;
        }
    }

    if let Some(variable_mods) = database
        .get_mut("variable_mods")
        .and_then(|v| v.as_object_mut())
    {
        for (_residue, v) in variable_mods.iter_mut() {
            if let Some(arr) = v.as_array_mut() {
                for entry in arr.iter_mut() {
                    resolve_mod_leaf(entry, &mut reverse)?;
                }
            }
        }
    }

    if !reverse.is_empty() {
        sage_core::unimod::set_reverse_table(reverse).map_err(|e| anyhow::anyhow!(e))?;
    }

    Ok(())
}

/// Same, for a full `Input`-shaped config value (`database` nested under
/// its own key) -- used by `Input::load`.
fn resolve_unimod_refs(value: &mut serde_json::Value) -> anyhow::Result<()> {
    if let Some(database) = value.get_mut("database") {
        resolve_unimod_refs_in_database(database)?;
    }
    Ok(())
}

impl Input {
    pub fn from_arguments(matches: ArgMatches) -> anyhow::Result<Self> {
        // Must happen before `Input::load` below: `--unimod-db-path`
        // overrides the embedded default table `static_mods`/
        // `variable_mods`'s `UNIMOD:<id>` references get resolved against.
        if let Some(unimod_db_path) = matches.get_one::<String>("unimod-db-path") {
            sage_core::unimod::set_active_table_from_path(std::path::Path::new(unimod_db_path))
                .map_err(|e| anyhow::anyhow!(e))?;
        }

        let path = matches
            .get_one::<String>("parameters")
            .expect("required parameters");
        let mut input = Input::load(path)
            .with_context(|| format!("Failed to read parameters from `{path}`"))?;

        // Handle JSON configuration overrides
        if let Some(output_directory) = matches.get_one::<String>("output_directory") {
            input.output_directory = Some(output_directory.into());
        }
        if let Some(fasta) = matches.get_one::<String>("fasta") {
            input.database.fasta = Some(fasta.into());
        }
        if let Some(mzml_paths) = matches.get_many::<String>("mzml_paths") {
            input.mzml_paths = Some(mzml_paths.into_iter().map(|p| p.into()).collect());
        }
        if let Some(pmsms) = matches.get_one::<String>("pmsms") {
            input.pmsms = Some(pmsms.into());
        }
        if let Some(precursors) = matches.get_one::<String>("precursors") {
            input.precursors = Some(precursors.into());
        }
        if let Some(predicted_properties) = matches.get_one::<String>("predicted-properties") {
            input.predicted_properties = Some(predicted_properties.into());
        }

        if let Some(write_pin) = matches.get_one::<bool>("write-pin").copied() {
            input.write_pin = Some(write_pin);
        }

        if let Some(write_report) = matches.get_one::<bool>("write-report").copied() {
            input.write_report = Some(write_report);
        }

        if let Some(annotate_matches) = matches.get_one::<bool>("annotate-matches").copied() {
            input.annotate_matches = Some(annotate_matches);
        }

        // avoid to later panic if these parameters are not set (but doesn't check if files exist)

        ensure!(
            input.database.fasta.is_some(),
            "`database.fasta` must be set. For more information try '--help'"
        );

        let pmsms_flags_given = [&input.pmsms, &input.precursors]
            .iter()
            .filter(|o| o.is_some())
            .count();
        ensure!(
            pmsms_flags_given == 0 || pmsms_flags_given == 2,
            "`--pmsms` and `--precursors` must be given together"
        );

        ensure!(
            pmsms_flags_given == 2
                || input
                    .mzml_paths
                    .as_ref()
                    .map(|p| p.len())
                    .unwrap_or_default()
                    > 0,
            "`mzml_paths` must be set, or `--pmsms`/`--precursors` given together. \
             For more information try '--help'"
        );

        Ok(input)
    }

    pub fn load<S: AsRef<str>>(path: S) -> anyhow::Result<Self> {
        let mut value: serde_json::Value =
            sage_cloudpath::util::read_json(path).map_err(anyhow::Error::from)?;
        resolve_unimod_refs(&mut value)?;
        Ok(serde_json::from_value(value)?)
    }

    fn check_mass_tolerances(tolerance: &Tolerance) {
        let (lo, hi) = match tolerance {
            Tolerance::Ppm(lo, hi) => (*lo, *hi),
            Tolerance::Pct(lo, hi) => {
                log::warn!(
                    "Pct tolerances are very rarely used for mass tolerances, did you mean ppm?"
                );
                (*lo, *hi)
            }
            Tolerance::Da(lo, hi) => (*lo, *hi),
        };
        if hi.abs() > lo.abs() {
            log::warn!(
                "Tolerances are applied to experimental masses, not theoretical: [{}, {}]",
                lo,
                hi
            );
        }
        if lo > 0.0 {
            log::warn!(
                "The `left` tolerance should probably be negative, for example: [{}, {}]",
                -lo,
                hi.abs()
            )
        }
        if hi < 0.0 {
            log::warn!(
                "The `right` tolerance should probably be positive, for example: [{}, {}]",
                -lo.abs(),
                hi
            )
        }
    }

    /// `rt_tol_sec`'s spline *values* are user-facing in seconds;
    /// `ProcessedSpectrum::scan_start_time` (the internal RT representation,
    /// see `spectrum.rs`) is in minutes. Only the values need conversion —
    /// the spline's own grid is already in `scan_start_time`'s native
    /// minutes (it's evaluated against observed RT directly, no unit
    /// mismatch there).
    fn rt_tol_sec_to_minutes(tolerance: ValueTolSpline) -> ValueTolSpline {
        let to_minutes = |mut spline: sage_core::spline::LinearSpline| {
            for v in spline.values.iter_mut() {
                *v /= 60.0;
            }
            spline
        };
        ValueTolSpline {
            lo: to_minutes(tolerance.lo),
            hi: to_minutes(tolerance.hi),
        }
    }

    pub fn build(mut self) -> anyhow::Result<Search> {
        let database = self.database.make_parameters();

        Self::check_mass_tolerances(&self.fragment_tol);
        Self::check_mass_tolerances(&self.precursor_tol);

        if self.predicted_properties.is_some()
            && (self.rt_tol_sec.is_none() || self.mobility_tol.is_none())
        {
            anyhow::bail!(
                "`predicted_properties` file supplied but `rt_tol_sec`/`mobility_tol` are not \
                 both configured — RT/IIM candidate filtering requires all three together. \
                 Either set both tolerances, or remove `predicted_properties`."
            );
        }
        if let Some(spline) = &self.rt_tol_sec {
            spline
                .validate()
                .map_err(|e| anyhow::anyhow!("invalid `rt_tol_sec`: {e}"))?;
        }
        if let Some(spline) = &self.mobility_tol {
            spline
                .validate()
                .map_err(|e| anyhow::anyhow!("invalid `mobility_tol`: {e}"))?;
        }

        if let Some(spline) = &self.fragment_tol_spline {
            spline
                .validate()
                .map_err(|e| anyhow::anyhow!("invalid `fragment_tol_spline`: {e}"))?;
            log::warn!(
                "Both `fragment_tol` and `fragment_tol_spline` are set — \
                 `fragment_tol_spline` takes over fragment matching entirely \
                 (including outside its grid range, via flat extrapolation); \
                 `fragment_tol` is only used for its own sanity-check warnings \
                 above and is otherwise unused."
            );
        }

        if let Some(isotope_errors) = self.isotope_errors {
            if isotope_errors.0 > isotope_errors.1 {
                log::error!("Minimum isotope_error value greater than maximum! Typical usage: `isotope_errors: [-1, 3]`");
                std::process::exit(1);
            }
        }
        if let Some(charges) = self.precursor_charge {
            if charges.0 > charges.1 {
                log::error!(
                    "Precursor charges should be specified [low, high], user provided: [{}, {}]",
                    charges.0,
                    charges.1
                );
                std::process::exit(1);
            }
        }

        if !self.predict_rt.unwrap_or(true)
            && self.quant.as_ref().and_then(|q| q.lfq).unwrap_or(false)
        {
            log::warn!(
                "`predict_rt: false` and `lfq: true` are incompatible. Setting `predict_rt: true`"
            );
            self.predict_rt = Some(true);
        }

        let pmsms_paths = match (&self.pmsms, &self.precursors) {
            (Some(pmsms), Some(precursors)) => Some(PmsmsPaths {
                pmsms: pmsms.into(),
                precursors: precursors.into(),
            }),
            _ => None,
        };

        // A pmsms pair still needs exactly one `mzml_paths` entry: the rest of the
        // runner keys per-file bookkeeping (file_id, output filenames) off that list.
        let mzml_paths = match &pmsms_paths {
            Some(p) => vec![sage_cloudpath::to_url(&p.pmsms.to_string_lossy())?],
            None => self
                .mzml_paths
                .expect("'mzml_paths' must be provided!")
                .iter()
                .map(|s| sage_cloudpath::to_url(s))
                .collect::<Result<Vec<_>, _>>()?,
        };

        let output_directory = match self.output_directory {
            Some(path) => {
                match Url::parse(&path) {
                    Ok(mut url) => {
                        // Valid URL, might still be a local directory that doesn't exist
                        if url.scheme() == "file" {
                            let path = url.to_file_path().expect("url scheme is file");
                            std::fs::create_dir_all(path)?;
                        }

                        if !url.path().ends_with("/") {
                            url.set_path(&format!("{}/", url.path()));
                        }
                        url
                    }
                    Err(_) => {
                        // Try to interpret as a local path
                        let path = std::path::Path::new(&path);
                        std::fs::create_dir_all(path)?;
                        Url::from_directory_path(path.canonicalize()?).expect("valid path")
                    }
                }
            }
            None => {
                let dir = std::env::current_dir()?;
                Url::from_directory_path(dir).expect("valid path")
            }
        };

        let score_type = self.score_type.unwrap_or(ScoreType::SageHyperScore);

        Ok(Search {
            version: clap::crate_version!().into(),
            database,
            quant: self.quant.map(Into::into).unwrap_or_default(),
            mzml_paths,
            pmsms_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            fragment_tol_spline: self.fragment_tol_spline,
            predicted_properties: self.predicted_properties,
            rt_tol: self.rt_tol_sec.map(Self::rt_tol_sec_to_minutes),
            mobility_tol: self.mobility_tol,
            report_psms: self.report_psms.unwrap_or(1),
            max_peaks: self.max_peaks.unwrap_or(150),
            min_peaks: self.min_peaks.unwrap_or(15),
            min_matched_peaks: self.min_matched_peaks.unwrap_or(4),
            max_fragment_charge: self.max_fragment_charge,
            annotate_matches: self.annotate_matches.unwrap_or(false),
            precursor_charge: self.precursor_charge.unwrap_or((2, 4)),
            override_precursor_charge: self.override_precursor_charge.unwrap_or(false),
            isotope_errors: self.isotope_errors.unwrap_or((0, 0)),
            deisotope: self.deisotope.unwrap_or(true),
            chimera: self.chimera.unwrap_or(false),
            wide_window: self.wide_window.unwrap_or(false),
            predict_rt: self.predict_rt.unwrap_or(true),
            output_paths: Vec::new(),
            write_pin: self.write_pin.unwrap_or(false),
            bruker_config: self.bruker_config.unwrap_or_default(),
            write_report: self.write_report.unwrap_or(false),
            protein_grouping: self.protein_grouping.unwrap_or(true),
            protein_grouping_peptide_fdr: self.protein_grouping_peptide_fdr.unwrap_or(0.01),
            score_type,
        })
    }
}

#[cfg(test)]
mod test {

    use super::{resolve_unimod_refs, Input};
    use sage_core::{database::EnzymeBuilder, enzyme::EnzymeParameters, mass::Tolerance};

    #[test]
    fn deserialize_enzyme_builder() -> Result<(), serde_json::Error> {
        let a: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
        }))?;
        let b: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "P",
        }))?;
        let c: EnzymeBuilder = serde_json::from_value(serde_json::json!({
            "cleave_at": "KR",
            "restrict": "",
        }))?;

        let a: EnzymeParameters = a.into();
        let b: EnzymeParameters = b.into();
        let c: EnzymeParameters = c.into();

        assert_eq!(a.enzyme.map(|e| e.skip_suffix), Some([false; 26]));
        {
            let mut expected = [false; 26];
            expected[(b'P' - b'A') as usize] = true;
            assert_eq!(b.enzyme.map(|e| e.skip_suffix), Some(expected));
        }
        assert_eq!(c.enzyme.map(|e| e.skip_suffix), Some([false; 26]));

        Ok(())
    }

    /// Minimal `log::Log` implementor that captures messages into a `Vec`,
    /// so tests can assert on a specific `log::warn!` without depending on
    /// a log-capturing crate. `log::set_logger` only succeeds once per
    /// process; safe here since no other test in this binary installs one.
    struct CapturingLogger {
        messages: std::sync::Mutex<Vec<String>>,
    }

    impl log::Log for CapturingLogger {
        fn enabled(&self, _metadata: &log::Metadata) -> bool {
            true
        }
        fn log(&self, record: &log::Record) {
            self.messages
                .lock()
                .unwrap()
                .push(record.args().to_string());
        }
        fn flush(&self) {}
    }

    fn install_capturing_logger() -> &'static CapturingLogger {
        static LOGGER: std::sync::OnceLock<CapturingLogger> = std::sync::OnceLock::new();
        let logger = LOGGER.get_or_init(|| CapturingLogger {
            messages: std::sync::Mutex::new(Vec::new()),
        });
        let _ = log::set_logger(logger);
        log::set_max_level(log::LevelFilter::Warn);
        logger
    }

    fn mk_input_json(fragment_tol_spline: Option<serde_json::Value>) -> serde_json::Value {
        // `Input::build` resolves `mzml_paths`/`database.fasta` to real
        // filesystem paths (canonicalize), so both must exist — reuse the
        // fixtures already committed for `sage-cli`'s other tests.
        let mut json = serde_json::json!({
            "database": {"fasta": "../../tests/Q99536.fasta"},
            "precursor_tol": {"ppm": [-10.0, 10.0]},
            "fragment_tol": {"ppm": [-10.0, 10.0]},
            "mzml_paths": ["../../tests/LQSRPAAPPAPGPGQLTLR.mzML"],
        });
        if let Some(spline) = fragment_tol_spline {
            json["fragment_tol_spline"] = spline;
        }
        json
    }

    fn flat_spline_json() -> serde_json::Value {
        serde_json::json!({
            "ppm_lo": {"grid_start": 0.0, "grid_step": 100.0, "values": [-10.0, -10.0]},
            "ppm_hi": {"grid_start": 0.0, "grid_step": 100.0, "values": [10.0, 10.0]},
        })
    }

    // Both scenarios live in one #[test] (rather than two) because they
    // share one process-global logger — cargo runs tests in the same
    // binary in parallel by default, and two tests independently
    // clear()-ing/reading that shared buffer would race.
    #[test]
    fn fragment_tol_spline_warning_only_fires_when_spline_is_set() {
        let logger = install_capturing_logger();

        logger.messages.lock().unwrap().clear();
        let input: Input = serde_json::from_value(mk_input_json(None)).unwrap();
        input.build().expect("valid config should build");
        assert!(
            !logger
                .messages
                .lock()
                .unwrap()
                .iter()
                .any(|m| m.contains("fragment_tol_spline")),
            "expected no fragment_tol_spline warning when spline is unset"
        );

        logger.messages.lock().unwrap().clear();
        let input: Input =
            serde_json::from_value(mk_input_json(Some(flat_spline_json()))).unwrap();
        input
            .build()
            .expect("both set is a warning, not a build error");
        let messages = logger.messages.lock().unwrap();
        assert!(
            messages
                .iter()
                .any(|m| m.contains("fragment_tol") && m.contains("fragment_tol_spline")),
            "expected a warning naming both `fragment_tol` and `fragment_tol_spline`, got: {:?}",
            *messages
        );
    }

    fn mk_predicted_properties_json(
        predicted_properties: Option<&str>,
        rt_tol_sec: Option<serde_json::Value>,
        mobility_tol: Option<serde_json::Value>,
    ) -> serde_json::Value {
        let mut json = mk_input_json(None);
        if let Some(path) = predicted_properties {
            json["predicted_properties"] = serde_json::json!(path);
        }
        if let Some(rt) = rt_tol_sec {
            json["rt_tol_sec"] = rt;
        }
        if let Some(im) = mobility_tol {
            json["mobility_tol"] = im;
        }
        json
    }

    /// A flat (value-independent) `ValueTolSpline` JSON blob, `[lo, hi]`
    /// everywhere -- 2 nodes with identical values, same shape used by the
    /// real `rt_tol_sec`/`mobility_tol` robust-flat-window fits.
    fn flat_value_tol_spline_json(lo: f64, hi: f64) -> serde_json::Value {
        serde_json::json!({
            "lo": {"grid_start": 0.0, "grid_step": 1.0, "values": [lo, lo]},
            "hi": {"grid_start": 0.0, "grid_step": 1.0, "values": [hi, hi]},
        })
    }

    #[test]
    fn predicted_properties_without_tolerances_errors() {
        let json = mk_predicted_properties_json(Some("predictions.parquet"), None, None);
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("missing tolerances should error"),
        };
        assert!(
            err.to_string().contains("rt_tol_sec") && err.to_string().contains("mobility_tol"),
            "expected error naming both `rt_tol_sec` and `mobility_tol`, got: {err}"
        );
    }

    #[test]
    fn predicted_properties_with_only_rt_tol_errors() {
        let json = mk_predicted_properties_json(
            Some("predictions.parquet"),
            Some(flat_value_tol_spline_json(-30.0, 30.0)),
            None,
        );
        let input: Input = serde_json::from_value(json).unwrap();
        assert!(
            input.build().is_err(),
            "partial predicted-properties config (rt_tol_sec set, mobility_tol unset) should error"
        );
    }

    #[test]
    fn predicted_properties_with_both_tolerances_converts_seconds_to_minutes() {
        let json = mk_predicted_properties_json(
            Some("predictions.parquet"),
            Some(flat_value_tol_spline_json(-30.0, 30.0)),
            Some(flat_value_tol_spline_json(-0.01, 0.01)),
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let search = input.build().expect("fully configured should build");
        assert_eq!(
            search.rt_tol.as_ref().unwrap().tolerance_at(0.0),
            Tolerance::Da(-0.5, 0.5),
            "rt_tol_sec values are in seconds, scan_start_time in minutes -- /60 at build time"
        );
        assert_eq!(
            search.mobility_tol.as_ref().unwrap().tolerance_at(0.0),
            Tolerance::Da(-0.01, 0.01),
            "mobility_tol is unitless (1/K0), no conversion"
        );
        assert_eq!(
            search.predicted_properties.as_deref(),
            Some("predictions.parquet")
        );
    }

    #[test]
    fn predicted_properties_with_invalid_rt_tol_spline_errors() {
        let mut bad_spline = flat_value_tol_spline_json(-30.0, 30.0);
        bad_spline["lo"]["grid_step"] = serde_json::json!(0.0);
        let json = mk_predicted_properties_json(
            Some("predictions.parquet"),
            Some(bad_spline),
            Some(flat_value_tol_spline_json(-0.01, 0.01)),
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("non-positive grid_step in rt_tol_sec should be rejected"),
        };
        assert!(
            err.to_string().contains("rt_tol_sec"),
            "expected error naming `rt_tol_sec`, got: {err}"
        );
    }

    #[test]
    fn no_predicted_properties_no_tolerance_error() {
        let json = mk_input_json(None);
        let input: Input = serde_json::from_value(json).unwrap();
        input
            .build()
            .expect("no predicted_properties, no tolerances required, should build fine");
    }

    /// Only test in this binary that resolves a `UNIMOD:<id>` reference --
    /// `sage_core::unimod`'s reverse table is a process-global `OnceLock`,
    /// settable once; a second test doing the same would race/conflict
    /// with cargo's default parallel-in-one-process test execution (same
    /// reason `fragment_tol_spline_warning_only_fires_when_spline_is_set`
    /// above combines two scenarios into one test).
    #[test]
    fn unimod_reference_resolves_to_embedded_table_mass() {
        let mut json = mk_input_json(None);
        json["database"]["static_mods"] = serde_json::json!({"C": "UNIMOD:4"});
        resolve_unimod_refs(&mut json).expect("UNIMOD:4 should resolve against the embedded table");

        let input: Input = serde_json::from_value(json).unwrap();
        let search = input.build().expect("valid config should build");
        assert!(
            search
                .database
                .static_mods
                .values()
                .any(|&m| (m - 57.021464).abs() < 1e-4),
            "expected static_mods to contain Carbamidomethyl's exact mass, got: {:?}",
            search.database.static_mods
        );
    }

    #[test]
    fn unknown_unimod_id_errors() {
        let mut json = mk_input_json(None);
        json["database"]["static_mods"] = serde_json::json!({"C": "UNIMOD:999999999"});
        let err = resolve_unimod_refs(&mut json).expect_err("unknown UNIMOD id should error");
        assert!(
            err.to_string().contains("999999999"),
            "expected error naming the missing id, got: {err}"
        );
    }

    #[test]
    fn non_unimod_string_mod_errors_clearly() {
        let mut json = mk_input_json(None);
        json["database"]["static_mods"] = serde_json::json!({"C": "not-a-unimod-ref"});
        let err = resolve_unimod_refs(&mut json).expect_err("non-UNIMOD string should error");
        assert!(
            err.to_string().contains("UNIMOD"),
            "expected error mentioning the `UNIMOD:<id>` requirement, got: {err}"
        );
    }

    #[test]
    fn plain_float_coinciding_with_known_mod_errors() {
        // Same value real job configs actually write for Carbamidomethyl
        // (see plans/) -- close to, but not bit-identical to, the
        // canonical 57.021464.
        let mut json = mk_input_json(None);
        json["database"]["static_mods"] = serde_json::json!({"C": 57.0216});
        let err = resolve_unimod_refs(&mut json)
            .expect_err("a float coinciding with a real Unimod entry should error");
        assert!(
            err.to_string().contains("UNIMOD:4"),
            "expected error naming UNIMOD:4 (Carbamidomethyl), got: {err}"
        );
    }

    #[test]
    fn plain_float_not_coinciding_with_anything_is_unaffected() {
        let mut json = mk_input_json(None);
        // Not close to any real Unimod entry.
        json["database"]["static_mods"] = serde_json::json!({"C": 12345.6789});
        resolve_unimod_refs(&mut json).expect("a genuinely novel mass should pass through");

        let input: Input = serde_json::from_value(json).unwrap();
        let search = input.build().expect("valid config should build");
        assert!(search
            .database
            .static_mods
            .values()
            .any(|&m| (m - 12345.6789).abs() < 1e-2));
    }
}
