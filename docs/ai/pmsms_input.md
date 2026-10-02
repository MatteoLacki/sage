# pmsms/precursors input

Core reader: `crates/sage-cloudpath/src/pmsms.rs`. Inputs, given as explicit
paths (never assumed to live together in one fixed-filename directory):

- `pmsms.mmappet/` — schema.txt + one `<idx>.bin` file per column (see
  `git/mmappet`'s Python writer), columns looked up by name. Needs
  `tof`(uint32) and `intensity`(uint32). An `mz` column, if present, is
  ignored with a warning.
- `tof2mz` table — an mmappet with one float32 column `mz` indexed by tof;
  fragment m/z is `table[tof]`. A tof outside the table is an error
  (`PmsmsError::TofOutOfRange`).
- precursors, either `.parquet` (via the `parquet` crate) or `.mmappet/`
  (schema.txt + one `<idx>.bin` file per column, columns looked up by name —
  see `git/mmappet`'s Python writer for the on-disk format). Needs 7 columns:
  `precursor_idx`(u64), `mz`(f64), `rt`(f64), `inv_ion_mobility`(f64),
  `charges`(i64), `fragment_spectrum_start`(u64), `fragment_event_cnt`(u64).
- optionally a fragment m/z recalibration, `--mz-recalibration` (below).

CLI: **`--pmsms <path> --precursors <path> --tof2mz <path>`**, all three
together — validated in `crates/sage-cli/src/input.rs::Input::from_arguments`.
Used instead of the positional spectra-paths argument; `PmsmsPaths`
(`crates/sage-cloudpath/src/util.rs`) carries the paths through `Search`. The
runner's `read_processed_spectra` branches on `parameters.pmsms_paths` to call
`pmsms::parse_processed`, bypassing the normal `FileFormat`-suffix dispatch.
`necromerge2`'s `run_sage` necroflow rule
(`git/ionmaidentools/src/ionmaidentools/pipelines.py`) passes exactly these.

There used to be two more inputs, both removed 2026-10 when the search path
stopped materializing fragment m/z (necromerge2
`plans/fragment_mz_correction_in_sage.md`): reading the pmsms' own `mz` column
when `--tof2mz` was absent, and a positional `<dir>.pmsms` directory with
`pmsms.mmappet/` and `precursors.parquet` at fixed names (`read_pmsms`, plus the
caller-less `read_pmsms_explicit`).

## `--tof2mz` table (2026-09, mandatory since 2026-10)

Why it was added: the recalibration search can run on raw MS2 spectra in place.
A precursor's `top_probable_frame_scan` spectrum is one raw MS2 spectrum, stored
contiguously and tof-sorted in `tdf2ms2`'s `events.ms2/data.mmappet` (tof,
intensity only), so a precursor table whose `fragment_spectrum_start`/
`fragment_event_cnt` point at that (frame, scan) slice, plus the dense tof -> m/z
table `materialize_pmsms_mz` itself uses, reproduces the mkpmsms-built search
exactly with no mkpmsms run and no m/z column (F9477: all 47,917 PSMs and 235,629
matched fragments identical; necromerge2
`git/pipeline_analysis/docs/ai/recalibration_paths.md`). Precursor slices may
overlap (several precursors sharing a spectrum). Since 2026-10 the final search
reads the pmsms the same way, with the recalibration applied in the reader (below).

Any pmsms read error (including the recalibration errors below) is only
logged: the runner turns it into an empty spectrum list, so SAGE searches 0
spectra and exits 0.

## Fragment m/z recalibration, `--mz-recalibration` (2026-10)

`--mz-recalibration <dir.mzcalib>`: timstofu's `MzRecalibration` artifact
(`bias.mmappet` plus one `<dim>/{grid,values}.mmappet` curve per dimension), here
with exactly one dimension, `mz`. Each fragment m/z read as `table[tof]` becomes

    mz / (1 + (bias + f_mz(mz) + fragment_shift_ppm) * 1e-6)

computed in f64 and rounded to f32, where `fragment_shift_ppm`(f64) is a column of
the precursors table: one constant per precursor. `FragmentMzRecalibration`
(`pmsms.rs`) loads and applies it per precursor in `parse_with`.

Errors: a `.mzcalib` with any dimension set other than `{mz}` (e.g. a leftover `rt`
curve, which would otherwise be silently skipped) is
`PmsmsError::MzRecalibrationDims`; a precursors table without
`fragment_shift_ppm` is `MissingColumn("fragment_shift_ppm")`. Without the flag
the column is ignored.

Why it is in the reader: the correction is tiny (2,000 grid points), while
materializing it upstream meant writing a 6.5 GB `mz` column (F9477) twice per run,
once for `table[tof]` and once recalibrated. Materialization now runs only for
exports (mzML/MGF/the SAGE-to-pmsms mapping), in timstofu's
`materialize_pmsms_mz --mz-recalibration ... --precursors ...`.

Why the per-precursor part is a column and not a curve SAGE evaluates: the fit is
`bias + f_mz(mz) + f_rt(precursor rt)`, and `f_rt` is constant per precursor.
searchops' `recalibrate_pmsms_mz` evaluates it once from the precursors' raw `rt`
and writes `fragment_shift_ppm = bias + f_rt(rt)` (the `.mzcalib` then carries
`bias = 0`). Evaluating `f_rt` in SAGE would have been wrong in necromerge2's
mode 3, where the RT correction overwrites the precursors' `rt` column with
corrected RT before the final search, while `f_rt` was fitted on raw RT. The
column also leaves the per-precursor model free to change without touching SAGE:
e.g. a 2-D `f(rt, 1/K0)`, or a free per-precursor intercept.

Bit-identity with the Python export path: `EvenlySpacedLinearSpline` mirrors
timstofu's `EvenlySpacedLinearSpline.njit_evaluator` operation for operation
(`spacing = (x_max - x_min) / (n - 1)`, `inv_spacing = 1 / spacing`, flat clamp,
truncating index), and `correct` mirrors `apply_fragment_mz_recalibration_blocks`'
`(bias + f_mz) + shift` order. SAGE's own `LinearSpline` (f32, `grid_step`) is not
reused: it would round differently. `recalibration_matches_timstofu_bit_for_bit`
pins both against values generated with the numba code.

## Borrowed fragment input and processed storage

The explicit-path CLI input calls `parse_processed`: each precursor's m/z
(looked up as `table[tof]` into a per-spectrum vector, then recalibrated in
place) and mapped `intensity: &[u32]` range feed `SpectrumProcessor` directly,
in parallel, while the reader owns the mappings. No intermediate per-spectrum
float-intensity vectors are allocated. Ordered parallel collection
preserves precursor order; returned processed spectra own their storage and
remain valid after the mappings close.

Intensity conversion to `f32` happens before every comparison and accumulation,
preserving the legacy rounding behavior even when distinct `u32` values round
to the same float. Deisotoping and top-N selection keep their existing ordering.
Final `PeakColumns` hold separate mass and intensity vectors, allocated from the
exact retained count; observed charges remain a parallel vector when present.
For max_peaks=800, each final column contains at most 800 entries. Temporary
selection/deisotoping storage still scales with the raw spectrum size.

`pmsms::parse` retains the owned `RawSpectrum` API (tests use it as the reference
for `parse_processed`); it copies raw columns.
Non-mobility processed spectra use `PeakColumns`; mobility spectra keep `Vec<IMPeak>`.
See [F9477 measurements](simd.md) for timing and memory validation.

## Per-precursor ppm tolerance

The precursors table accepts two further **optional** columns: `ppm_tol_lo`,
`ppm_tol_hi` (f64, signed — lo negative, hi positive). When both are present
for a row, that precursor is searched with `Tolerance::Ppm(ppm_tol_lo,
ppm_tol_hi)` instead of the run's global `--precursor_tol`; when either is
missing (the common case — no upstream writer populates them yet), behavior
is unchanged from before this existed. Design/history:
`necromerge2`'s `plans/per_precursor_ppm_tolerance.md`.

This rides on `Precursor::isolation_window: Option<Tolerance>`
(`crates/sage/src/spectrum.rs`), which pre-existing wide-window/DIA search
already used per-spectrum; `Precursor::effective_precursor_tol` is the
`unwrap_or`-onto-global-tol helper, used by both the wide-window and
standard search branches in `crates/sage/src/scoring.rs::initial_hits`.
Column reading lives in `crates/sage-cloudpath/src/pmsms.rs` (parquet and
mmappet paths both use an *optional* column lookup — missing column is
`None`, not an error, unlike the 7 required columns).

## Per-fragment ppm tolerance: mass-dependent spline tried, then removed (2026-09)

A `Scorer::fragment_tol_spline: Option<FragmentTolSpline>` existed briefly
(`crates/sage/src/spline.rs`, design in `necromerge2`'s
`plans/per_fragment_ppm_tolerance.md`), letting the fragment ppm window vary
with fragment mass instead of one flat `fragment_tol`. Dropped — never wired
into any real job config here — after real-data analysis (real F9477
confident-PSM fragment residuals, `git/searchops`'s `recalibrate_pmsms_mz`
fit) found two problems with it: it was evaluated at `peak.mass` (already
charge-corrected neutral mass) while the calibration it would have needed to
follow is fit on raw experimental *m/z*, a real mismatch; and the bigger
finding — the *existing* flat tolerance's own calibration was itself biased,
because the initial SAGE pass's narrow capture window (±15 ppm) truncates
the true residual distribution before `_select_tolerance` ever sees it
(confirmed: widening to ±25 ppm recovered ~4.4% more confident fragments and
revealed materially wider true tails, ~19–24 ppm vs the ±15 ppm window's own
edge). Fixing that truncation is the more consequential change; a
mass-dependent window was solving a smaller, secondary effect. `page_search`
still takes `fragment_tol` as an explicit per-call argument (general,
independently useful — e.g. `wide_window` mode), just no longer varies it
by mass.

Separately, even taken on its own merits, mass-dependence wasn't worth its
cost: per-m/z-bin residual spread (real F9477 confident targets) only
diverged meaningfully from the ~5–6 ppm core below ~300 m/z (~7–8 ppm
there, traced to real tof-quantization physics — a single tof step is
already ~10-16 ppm of quantization noise below 300 m/z, vs ~5-7 ppm through
the 500–1000 m/z bulk). That region is a real but minority slice of the
fragment population (~8% of confident fragments below 300 m/z in this
dataset) — not enough weight to justify a per-peak spline evaluation on
every fragment match, an extra config surface, and the correctness risk
that bit us above, for the whole run. (The high-m/z tail's apparent
widening, by contrast, wasn't well-supported at all — driven by bins down
to a few dozen points, noise rather than signal, not a reason to act on
either.)

## MS1 precursor intensity (optional, apex only) and MS2 intensity, 2026-09

The precursors table accepts one further **optional** column: `intensity`
(i64 in both the parquet and mmappet fixture layouts — not `uint32`, despite
the fragment-level `pmsms.mmappet` column of the same name being genuinely
`uint32`; the two live in different files and have never had matching
dtypes). When present, it becomes `Precursor::intensity` →
`Feature::ms1_intensity` → the PIN's `ln(ms1_intensity)` column; when absent,
`ms1_intensity` reports `0.0`, same optional-column convention as
`ppm_tol_lo`/`ppm_tol_hi` above.

This is **timstofu's apex-point MS1 reading** at the precursor's selected
(tof,urt,scan) position, not an isotope-envelope sum — `matched_intensity_pct`
is a real sum for MS2 because deisotoping folds each isotope satellite's
intensity into its monoisotopic root before scoring, but no equivalent
folding happens on the MS1 side before this value is captured. It reaches the
precursors table only because `git/mgfops`'s `cut_and_index_precursors.py`
merges the full precursor-cluster row (every column timstofu produced,
including its own `intensity` and `log10_intensity`) against the fragment
index — not a deliberate export. A true envelope-summed MS1 intensity would
need timstofu's `intensity_argmax` candidate-selection method (its
`box_intensity` output) wired through instead; the production config
(`correlation` method) never computes one. See `necromerge2`'s
`plans/isotope_envelope_scoring.md` for the deferred IsoSpec-based
envelope-comparison follow-on this is a prerequisite for.

`Feature::ms2_intensity` (`summed_b + summed_y`, already isotope-summed via
deisotoping — see above) existed on `Feature` before this work and needed no
reader change; it was simply never written to the PIN. Now emitted as
`ln(ms2_intensity)`, placed next to the other `ms2_*` PIN columns rather than
next to `ln(ms1_intensity)`, which sits with the other raw precursor-signal
columns (`retentiontime`/`ion_mobility`) instead — grouped by what each
column describes, not by when it was added.
