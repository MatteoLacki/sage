use anyhow::{ensure, Context};
use clap::ArgMatches;
use sage_cloudpath::tdf::BrukerProcessingConfig;
use sage_cloudpath::util::PmsmsPaths;
use sage_cloudpath::Url;
use sage_core::scoring::{RankingScore, ScoreType};
use sage_core::spline::{LinearSpline, ValueTolSpline};
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
    /// Predicted-RT/IIM candidate filtering — see `predicted_rt`/
    /// `predicted_iim` docs on [`Input`]. Independent of `predicted_iim`.
    /// `rt_tol`'s spline *values* here are already converted to minutes;
    /// its grid stays in `scan_start_time`'s native minutes too.
    pub predicted_rt: Option<String>,
    /// Independent of `predicted_rt` — see `predicted_iim` doc on [`Input`].
    pub predicted_iim: Option<String>,
    pub rt_tol: Option<ValueTolSpline>,
    pub mobility_tol: Option<ValueTolSpline>,
    /// Robust (MAD-based) scale of the `predicted_rt` residual as a
    /// function of observed `scan_start_time`, for normalizing
    /// `Feature::delta_rt_z2_external` into a z² LDA feature — see
    /// `plans/lda_external_rt_iim_features.md`. Already converted to
    /// minutes at build time, same as `rt_tol`'s spline values.
    pub rt_sigma: Option<LinearSpline>,
    /// Same as `rt_sigma`, for `predicted_iim` — unitless (1/K0), no
    /// conversion needed, same as `mobility_tol`.
    pub iim_sigma: Option<f32>,
    /// Path to the job-scoped `(sequence, charge, start, end)` pointer
    /// parquet (`git/featureprediction`'s `export_fragment_intensity_for_sage`).
    /// Feature-only (no hard filtering) — see `predicted_fragment_intensity_cache`
    /// and `docs/ai/predicted_fragment_intensity.md`. Must be set together
    /// with `predicted_fragment_intensity_cache`, or neither.
    pub predicted_fragment_intensity_index: Option<String>,
    /// Path to the shared `arrays.mmappet` directory (a subdirectory of
    /// `git/featureprediction`'s fragment-intensity `PredictionCache` —
    /// this fork never reads that cache's `index.sqlite3`/`write.lock`,
    /// only its arrays) that `predicted_fragment_intensity_index`'s
    /// `(start, end)` ranges address.
    pub predicted_fragment_intensity_cache: Option<String>,
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
    pub ranking_score: RankingScore,
}

#[derive(Deserialize)]
/// Input search parameters deserialized from JSON file
pub struct Input {
    pub database: Builder,
    pub precursor_tol: Tolerance,
    pub fragment_tol: Tolerance,
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
    /// Path to a parquet file of externally-predicted peptide RT (columns:
    /// sequence, rt), used to reject candidates whose predicted RT falls
    /// outside `rt_tol_sec` of the observed spectrum's RT. Requires
    /// `rt_tol_sec` to be set. Independent of `predicted_iim` — either,
    /// both, or neither may be given.
    pub predicted_rt: Option<String>,
    /// Path to a parquet file of externally-predicted peptide IIM (columns:
    /// sequence, charge, iim), used to reject candidates whose predicted
    /// IIM falls outside `mobility_tol` of the observed spectrum's IIM.
    /// Requires `mobility_tol` to be set. Independent of `predicted_rt`.
    pub predicted_iim: Option<String>,
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
    /// Robust (MAD-based) scale of the `predicted_rt` residual as a
    /// function of observed `scan_start_time` (seconds, same unit and grid
    /// convention as `rt_tol_sec` — converted to minutes at build time). A
    /// flat (value-independent) scale is just a 1- or 2-node spline with
    /// identical values, same convention as `rt_tol_sec`'s own flat case —
    /// there is no separate scalar type. Required alongside
    /// `predicted_rt`/`rt_tol_sec` — feeds `Feature::delta_rt_z2_external`,
    /// evaluated at each candidate's own `scan_start_time` rather than a
    /// single global scale, see `plans/lda_external_rt_iim_features.md` and
    /// `plans/rt_heteroscedastic_tolerance_spline.md`.
    pub rt_sigma_sec: Option<LinearSpline>,
    /// Same as `rt_sigma_sec`, for `predicted_iim`/`mobility_tol` — unitless
    /// (1/K0), no conversion needed.
    pub iim_sigma: Option<f32>,
    /// Path to a job-scoped `(sequence, charge, start, end)` pointer parquet
    /// (`git/featureprediction`'s `export_fragment_intensity_for_sage`),
    /// used purely to compute `Feature::ms2_entropy_similarity` -- no hard
    /// eviction, unlike `predicted_rt`/`predicted_iim` above. Must be given
    /// together with `predicted_fragment_intensity_cache`, or neither.
    pub predicted_fragment_intensity_index: Option<String>,
    /// Path to the shared `arrays.mmappet` directory
    /// `predicted_fragment_intensity_index`'s `(start, end)` ranges address.
    pub predicted_fragment_intensity_cache: Option<String>,
    pub bruker_config: Option<BrukerProcessingConfig>,
    pub protein_grouping: Option<bool>,
    pub protein_grouping_peptide_fdr: Option<f32>,

    pub annotate_matches: Option<bool>,
    pub write_pin: Option<bool>,
    pub write_report: Option<bool>,
    pub score_type: Option<ScoreType>,
    pub ranking_score: Option<RankingScore>,
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
        if let Some(predicted_rt) = matches.get_one::<String>("predicted-rt") {
            input.predicted_rt = Some(predicted_rt.into());
        }
        if let Some(predicted_iim) = matches.get_one::<String>("predicted-iim") {
            input.predicted_iim = Some(predicted_iim.into());
        }
        if let Some(path) = matches.get_one::<String>("predicted-fragment-intensity-index") {
            input.predicted_fragment_intensity_index = Some(path.into());
        }
        if let Some(path) = matches.get_one::<String>("predicted-fragment-intensity-cache") {
            input.predicted_fragment_intensity_cache = Some(path.into());
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

    /// `rt_tol_sec`/`rt_sigma_sec`'s spline *values* are user-facing in
    /// seconds; `ProcessedSpectrum::scan_start_time` (the internal RT
    /// representation, see `spectrum.rs`) is in minutes. Only the values
    /// need conversion — the spline's own grid is already in
    /// `scan_start_time`'s native minutes (it's evaluated against observed
    /// RT directly, no unit mismatch there). Generic over spline length, so
    /// this applies unchanged to `rt_sigma_sec` now that it's spline-shaped
    /// too, not just `rt_tol_sec`'s two nodes.
    fn spline_secs_to_minutes(mut spline: LinearSpline) -> LinearSpline {
        for v in spline.values.iter_mut() {
            *v /= 60.0;
        }
        spline
    }

    fn rt_tol_sec_to_minutes(tolerance: ValueTolSpline) -> ValueTolSpline {
        ValueTolSpline {
            lo: Self::spline_secs_to_minutes(tolerance.lo),
            hi: Self::spline_secs_to_minutes(tolerance.hi),
        }
    }

    pub fn build(mut self) -> anyhow::Result<Search> {
        let database = self.database.make_parameters();

        Self::check_mass_tolerances(&self.fragment_tol);
        Self::check_mass_tolerances(&self.precursor_tol);

        if !(self.predicted_rt.is_some() == self.rt_tol_sec.is_some()
            && self.rt_tol_sec.is_some() == self.rt_sigma_sec.is_some())
        {
            anyhow::bail!(
                "`predicted_rt`, `rt_tol_sec`, and `rt_sigma_sec` must all be set together \
                 (or all omitted) — either set all three, or remove whichever are set."
            );
        }
        if !(self.predicted_iim.is_some() == self.mobility_tol.is_some()
            && self.mobility_tol.is_some() == self.iim_sigma.is_some())
        {
            anyhow::bail!(
                "`predicted_iim`, `mobility_tol`, and `iim_sigma` must all be set together \
                 (or all omitted) — either set all three, or remove whichever are set."
            );
        }
        if self.predicted_fragment_intensity_index.is_some()
            != self.predicted_fragment_intensity_cache.is_some()
        {
            anyhow::bail!(
                "`predicted_fragment_intensity_index` and `predicted_fragment_intensity_cache` \
                 must both be set together (or both omitted)."
            );
        }
        if let Some(spline) = &self.rt_tol_sec {
            spline
                .validate()
                .map_err(|e| anyhow::anyhow!("invalid `rt_tol_sec`: {e}"))?;
        }
        if let Some(spline) = &self.rt_sigma_sec {
            spline
                .validate()
                .map_err(|e| anyhow::anyhow!("invalid `rt_sigma_sec`: {e}"))?;
        }
        if let Some(spline) = &self.mobility_tol {
            spline
                .validate()
                .map_err(|e| anyhow::anyhow!("invalid `mobility_tol`: {e}"))?;
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
        let ranking_score = self.ranking_score.unwrap_or(RankingScore::CombinedScore);

        Ok(Search {
            version: clap::crate_version!().into(),
            database,
            quant: self.quant.map(Into::into).unwrap_or_default(),
            mzml_paths,
            pmsms_paths,
            output_directory,
            precursor_tol: self.precursor_tol,
            fragment_tol: self.fragment_tol,
            predicted_rt: self.predicted_rt,
            predicted_iim: self.predicted_iim,
            rt_tol: self.rt_tol_sec.map(Self::rt_tol_sec_to_minutes),
            mobility_tol: self.mobility_tol,
            rt_sigma: self.rt_sigma_sec.map(Self::spline_secs_to_minutes),
            iim_sigma: self.iim_sigma,
            predicted_fragment_intensity_index: self.predicted_fragment_intensity_index,
            predicted_fragment_intensity_cache: self.predicted_fragment_intensity_cache,
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
            ranking_score,
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

    fn mk_input_json() -> serde_json::Value {
        // `Input::build` resolves `mzml_paths`/`database.fasta` to real
        // filesystem paths (canonicalize), so both must exist — reuse the
        // fixtures already committed for `sage-cli`'s other tests.
        serde_json::json!({
            "database": {"fasta": "../../tests/Q99536.fasta"},
            "precursor_tol": {"ppm": [-10.0, 10.0]},
            "fragment_tol": {"ppm": [-10.0, 10.0]},
            "mzml_paths": ["../../tests/LQSRPAAPPAPGPGQLTLR.mzML"],
        })
    }

    fn mk_predicted_json(
        predicted_rt: Option<&str>,
        rt_tol_sec: Option<serde_json::Value>,
        predicted_iim: Option<&str>,
        mobility_tol: Option<serde_json::Value>,
    ) -> serde_json::Value {
        mk_predicted_json_with_sigma(predicted_rt, rt_tol_sec, None, predicted_iim, mobility_tol, None)
    }

    /// Same as `mk_predicted_json`, plus `rt_sigma_sec`/`iim_sigma` — the
    /// two are now required alongside `rt_tol_sec`/`mobility_tol`
    /// respectively (`plans/lda_external_rt_iim_features.md`), so tests
    /// exercising the fully-configured/converts-successfully paths need a
    /// way to set them; tests exercising the "missing" validation errors
    /// use `mk_predicted_json` (which leaves them `None`) instead.
    /// `rt_sigma_sec` takes a JSON `LinearSpline` blob (see
    /// `flat_linear_spline_json`), not a bare scalar, now that it's
    /// RT-dependent (`plans/rt_heteroscedastic_tolerance_spline.md`);
    /// `iim_sigma` stays a bare scalar -- IIM is out of scope for that
    /// change.
    fn mk_predicted_json_with_sigma(
        predicted_rt: Option<&str>,
        rt_tol_sec: Option<serde_json::Value>,
        rt_sigma_sec: Option<serde_json::Value>,
        predicted_iim: Option<&str>,
        mobility_tol: Option<serde_json::Value>,
        iim_sigma: Option<f32>,
    ) -> serde_json::Value {
        let mut json = mk_input_json();
        if let Some(path) = predicted_rt {
            json["predicted_rt"] = serde_json::json!(path);
        }
        if let Some(rt) = rt_tol_sec {
            json["rt_tol_sec"] = rt;
        }
        if let Some(sigma) = rt_sigma_sec {
            json["rt_sigma_sec"] = sigma;
        }
        if let Some(path) = predicted_iim {
            json["predicted_iim"] = serde_json::json!(path);
        }
        if let Some(im) = mobility_tol {
            json["mobility_tol"] = im;
        }
        if let Some(sigma) = iim_sigma {
            json["iim_sigma"] = serde_json::json!(sigma);
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

    /// A flat (value-independent) `LinearSpline` JSON blob -- single-value
    /// grid, constant everywhere (`LinearSpline::eval`'s single-value case).
    /// Same role for `rt_sigma_sec` as `flat_value_tol_spline_json` has for
    /// `rt_tol_sec`, now that both are spline-shaped
    /// (`plans/rt_heteroscedastic_tolerance_spline.md`).
    fn flat_linear_spline_json(value: f64) -> serde_json::Value {
        serde_json::json!({"grid_start": 0.0, "grid_step": 1.0, "values": [value]})
    }

    #[test]
    fn predicted_rt_without_tolerance_errors() {
        let json = mk_predicted_json(Some("predicted_rt.parquet"), None, None, None);
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("predicted_rt without rt_tol_sec should error"),
        };
        assert!(
            err.to_string().contains("predicted_rt") && err.to_string().contains("rt_tol_sec"),
            "expected error naming `predicted_rt` and `rt_tol_sec`, got: {err}"
        );
    }

    #[test]
    fn predicted_iim_without_tolerance_errors() {
        let json = mk_predicted_json(None, None, Some("predicted_iim.parquet"), None);
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("predicted_iim without mobility_tol should error"),
        };
        assert!(
            err.to_string().contains("predicted_iim") && err.to_string().contains("mobility_tol"),
            "expected error naming `predicted_iim` and `mobility_tol`, got: {err}"
        );
    }

    #[test]
    fn predicted_rt_only_builds_successfully() {
        // RT-only filtering is a supported, independent mode -- no mobility_tol
        // needed at all. See plans/rt_iim_independent_dimensions.md.
        let json = mk_predicted_json_with_sigma(
            Some("predicted_rt.parquet"),
            Some(flat_value_tol_spline_json(-30.0, 30.0)),
            Some(flat_linear_spline_json(5.0)),
            None,
            None,
            None,
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let search = input
            .build()
            .expect("RT-only predicted-property filtering should be valid");
        assert_eq!(
            search.predicted_rt.as_deref(),
            Some("predicted_rt.parquet")
        );
        assert!(search.predicted_iim.is_none());
        assert!(search.mobility_tol.is_none());
        assert!(search.iim_sigma.is_none());
    }

    #[test]
    fn predicted_iim_only_builds_successfully() {
        let json = mk_predicted_json_with_sigma(
            None,
            None,
            None,
            Some("predicted_iim.parquet"),
            Some(flat_value_tol_spline_json(-0.01, 0.01)),
            Some(0.005),
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let search = input
            .build()
            .expect("IIM-only predicted-property filtering should be valid");
        assert_eq!(
            search.predicted_iim.as_deref(),
            Some("predicted_iim.parquet")
        );
        assert!(search.predicted_rt.is_none());
        assert!(search.rt_tol.is_none());
        assert!(search.rt_sigma.is_none());
    }

    #[test]
    fn predicted_rt_and_iim_together_converts_seconds_to_minutes() {
        let json = mk_predicted_json_with_sigma(
            Some("predicted_rt.parquet"),
            Some(flat_value_tol_spline_json(-30.0, 30.0)),
            Some(flat_linear_spline_json(6.0)),
            Some("predicted_iim.parquet"),
            Some(flat_value_tol_spline_json(-0.01, 0.01)),
            Some(0.005),
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
            search.rt_sigma.as_ref().unwrap().eval(0.0),
            0.1,
            "rt_sigma_sec is in seconds, same /60 conversion as rt_tol_sec's values"
        );
        assert_eq!(
            search.iim_sigma,
            Some(0.005),
            "iim_sigma is unitless (1/K0), no conversion"
        );
        assert_eq!(
            search.predicted_rt.as_deref(),
            Some("predicted_rt.parquet")
        );
        assert_eq!(
            search.predicted_iim.as_deref(),
            Some("predicted_iim.parquet")
        );
    }

    #[test]
    fn predicted_rt_with_invalid_rt_tol_spline_errors() {
        let mut bad_spline = flat_value_tol_spline_json(-30.0, 30.0);
        bad_spline["lo"]["grid_step"] = serde_json::json!(0.0);
        let json = mk_predicted_json_with_sigma(
            Some("predicted_rt.parquet"),
            Some(bad_spline),
            Some(flat_linear_spline_json(5.0)),
            None,
            None,
            None,
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
    fn predicted_rt_with_invalid_rt_sigma_spline_errors() {
        // rt_sigma_sec is now spline-shaped too (LinearSpline, not a bare
        // scalar) -- plans/rt_heteroscedastic_tolerance_spline.md -- so it
        // needs its own validate() check, mirroring rt_tol_sec's.
        let mut bad_sigma = flat_linear_spline_json(5.0);
        bad_sigma["grid_step"] = serde_json::json!(0.0);
        let json = mk_predicted_json_with_sigma(
            Some("predicted_rt.parquet"),
            Some(flat_value_tol_spline_json(-30.0, 30.0)),
            Some(bad_sigma),
            None,
            None,
            None,
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("non-positive grid_step in rt_sigma_sec should be rejected"),
        };
        assert!(
            err.to_string().contains("rt_sigma_sec"),
            "expected error naming `rt_sigma_sec`, got: {err}"
        );
    }

    #[test]
    fn predicted_rt_without_sigma_errors() {
        // rt_tol_sec set, rt_sigma_sec omitted -- the new three-way pairing
        // requirement (plans/lda_external_rt_iim_features.md) should reject
        // this the same way it already rejects predicted_rt without rt_tol_sec.
        let json = mk_predicted_json(
            Some("predicted_rt.parquet"),
            Some(flat_value_tol_spline_json(-30.0, 30.0)),
            None,
            None,
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("rt_tol_sec without rt_sigma_sec should error"),
        };
        assert!(
            err.to_string().contains("rt_sigma_sec"),
            "expected error naming `rt_sigma_sec`, got: {err}"
        );
    }

    #[test]
    fn predicted_iim_without_sigma_errors() {
        let json = mk_predicted_json(
            None,
            None,
            Some("predicted_iim.parquet"),
            Some(flat_value_tol_spline_json(-0.01, 0.01)),
        );
        let input: Input = serde_json::from_value(json).unwrap();
        let err = match input.build() {
            Err(e) => e,
            Ok(_) => panic!("mobility_tol without iim_sigma should error"),
        };
        assert!(
            err.to_string().contains("iim_sigma"),
            "expected error naming `iim_sigma`, got: {err}"
        );
    }

    #[test]
    fn no_predicted_rt_iim_no_tolerance_error() {
        let json = mk_input_json();
        let input: Input = serde_json::from_value(json).unwrap();
        input
            .build()
            .expect("no predicted_rt/predicted_iim, no tolerances required, should build fine");
    }

    /// Only test in this binary that resolves a `UNIMOD:<id>` reference --
    /// `sage_core::unimod`'s reverse table is a process-global `OnceLock`,
    /// settable once; a second test doing the same would race/conflict
    /// with cargo's default parallel-in-one-process test execution.
    #[test]
    fn unimod_reference_resolves_to_embedded_table_mass() {
        let mut json = mk_input_json();
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
        let mut json = mk_input_json();
        json["database"]["static_mods"] = serde_json::json!({"C": "UNIMOD:999999999"});
        let err = resolve_unimod_refs(&mut json).expect_err("unknown UNIMOD id should error");
        assert!(
            err.to_string().contains("999999999"),
            "expected error naming the missing id, got: {err}"
        );
    }

    #[test]
    fn non_unimod_string_mod_errors_clearly() {
        let mut json = mk_input_json();
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
        let mut json = mk_input_json();
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
        let mut json = mk_input_json();
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
