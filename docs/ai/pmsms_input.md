# pmsms/precursors input

Core reader: `crates/sage-cloudpath/src/pmsms.rs`. Two inputs, given as
explicit paths (never assumed to live together in one fixed-filename
directory):

- `pmsms.mmappet/` — schema.txt + one `<idx>.bin` file per column (see
  `git/mmappet`'s Python writer), columns looked up by name. Needs
  `mz`(float32) and `intensity`(uint32) — `mz` is materialized upstream by
  `timstofu`'s `materialize_pmsms_mz`/`recalibrate_pmsms_mz` (see necromerge2
  `plans/mzs_instead_tofs.md`); this crate never reads `tof` or looks up a
  separate tof2mz table.
- precursors, either `.parquet` (via the `parquet` crate) or `.mmappet/`
  (schema.txt + one `<idx>.bin` file per column, columns looked up by name —
  see `git/mmappet`'s Python writer for the on-disk format). Needs 7 columns:
  `precursor_idx`(u64), `mz`(f64), `rt`(f64), `inv_ion_mobility`(f64),
  `charges`(i64), `fragment_spectrum_start`(u64), `fragment_event_cnt`(u64).

Two ways to pass these two paths to the `sage` CLI:

1. **`--pmsms <path> --precursors <path>`** (preferred). Must be given
   together — validated in `crates/sage-cli/src/input.rs::Input::from_arguments`.
   Used instead of the positional spectra-paths argument; `PmsmsPaths`
   (`crates/sage-cloudpath/src/util.rs`) carries the two paths through
   `Search`, and the runner's per-file loop
   (`crates/sage-cli/src/runner.rs::process_chunk`) branches on
   `parameters.pmsms_paths` to call `read_pmsms_explicit` directly, bypassing
   the normal `FileFormat`-suffix dispatch.
2. **Positional `<dir>.pmsms`** (legacy, still works, unchanged behavior) — a
   directory whose name ends in `.pmsms`, containing `pmsms.mmappet/` and
   `precursors.parquet` at fixed names. Detected by
   `crates/sage-cloudpath/src/util.rs::FileFormat::from` and read via
   `read_pmsms`, which just joins the two fixed names and calls the same
   `pmsms::parse`.

`necromerge2`'s `run_sage` Snakemake-via-necroflow rule
(`git/ionmaidentools/src/ionmaidentools/pipelines.py`) uses option 1 — no
staging directory, no symlinks, no mmappet→parquet conversion needed anymore.

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
