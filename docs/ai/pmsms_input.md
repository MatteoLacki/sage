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

## Per-fragment ppm tolerance (mass-dependent, spline-based)

`Scorer::fragment_tol_spline: Option<FragmentTolSpline>`
(`crates/sage/src/spline.rs`) lets the fragment ppm window vary with
fragment mass instead of being one flat `fragment_tol` for the whole run —
given as two independent `LinearSpline`s (`ppm_lo`, `ppm_hi`), each a
piecewise-linear function sampled on an equally spaced grid. Extrapolation
outside the grid range is controlled by `LinearSpline::extrapolation`
(`Extrapolation::Flat` — the default, and `FragmentTolSpline`'s own current
usage — clamps to the nearest edge value; `Extrapolation::Linear` continues
the boundary segment's slope instead, added for `git/featureprediction`'s
Python `LinearSpline` port — see that repo's `AI.md` — which needed the
opposite default). Missing `extrapolation` in a JSON config deserializes to
`Flat`, so existing configs are unaffected. `None` fragment_tol_spline (the
default — no JSON config key at all) behaves exactly as before this existed.
Design/history: `necromerge2`'s `plans/per_fragment_ppm_tolerance.md`.

Configured via `Input`/`Search`'s `fragment_tol_spline` JSON field
(`crates/sage-cli/src/input.rs`), validated once in `Input::build()`
(non-empty grid, positive `grid_step`).

Evaluated **per observed peak**, at the peak's own mass, in
`crates/sage/src/scoring.rs`, at all three places fragment matching happens
against `self.fragment_tol`:
- `Scorer::matched_peaks_with_isotope` (the preliminary candidate search —
  observed peak known, spline evaluated at `peak.mass`),
- `Scorer::score_candidate` (the final per-candidate hyperscore pass, which
  walks *theoretical* ions and searches for a matching peak — spline
  evaluated at the theoretical `mz` instead, since the observed peak isn't
  known yet; a genuine match differs from theoretical by less than the
  tolerance width itself, so this is equivalent in practice),
- `Scorer::remove_matched_peaks` (chimera-mode peak stripping between
  passes).

**All three must stay in sync** — the preliminary pass alone finding a
candidate isn't sufficient; `score_candidate`'s independent re-match (using
a stale flat `self.fragment_tol`) will silently zero out `matched_b`/`matched_y`
and drop the candidate even though the preliminary pass found it. (This bit
us during implementation — the reachability integration test caught it.)

`IndexedQuery::page_search` (`crates/sage/src/database.rs`) takes the
tolerance as an explicit argument rather than reading a fixed one baked into
the query, specifically so callers can vary it per call.
