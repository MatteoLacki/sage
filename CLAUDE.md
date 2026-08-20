# CLAUDE.md — sage fork conventions

This is Matteo's patched fork of [lazear/sage](https://github.com/lazear/sage)
(`git@github.com:MatteoLacki/sage.git`, branch `master`), vendored into the
necromerge2 monorepo at `software/sage/devel_fixed`. It's an independent git
repository — `necromerge2`'s root repo has `software` in its `.gitignore`, and
this fork is not a submodule. Built via `cargo build --release --bin sage`;
`necromerge2`'s `Makefile` symlinks `software/sage/devel_fixed/sage →
target/release/sage`. The release profile uses `lto = "fat"`, so a full
release build takes ~1-2 minutes.

## Patches over upstream

- Fragment-reporting fix (see Makefile comment in the parent repo: "SAGE from
  github fork with fragment reporting fix").
- Native support for reading ionmaiden's `pmsms`/precursors binary formats
  directly, instead of only mzML/MGF/(Bruker) TDF.

## pmsms/precursors input

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

## Predicted RT/IIM candidate filtering

`--predicted-properties <path.parquet>` (columns `sequence, charge, rt, iim`)
loads a read-only `HashMap<(String, u8), (f32, f32)>` once at startup
(`Runner::predicted_properties`, `crates/sage-cli/src/runner.rs`), shared
across the parallel search same as `database`. Requires config fields
`rt_tol_sec` and `mobility_tol` to both be set; `Input::build()` rejects
`--predicted-properties` given without both. `Scorer::evict_rt_iim_mismatches`
(`crates/sage/src/scoring.rs`) evicts fragment-matched candidates whose
predicted RT/IIM falls outside the observed spectrum's window, right before
`trim_hits`. A candidate with no `(sequence, charge)` entry in the map is
left alone (permissive — avoids rejecting due to prediction-coverage gaps).
Design/history: `necromerge2`'s `plans/better_sage_filtering.md`.

Both fields are `ValueTolSpline` (`crates/sage/src/spline.rs`) — the same
two-independent-`LinearSpline` (`lo`/`hi`) shape as `FragmentTolSpline`,
evaluated against one observed value (`ProcessedSpectrum::scan_start_time`
for RT, `Precursor::inverse_ion_mobility` for IIM) rather than a flat
`Tolerance::Da`. Deliberately generic (not `RtTolSpline`/`MobilityTolSpline`)
since RT and IIM tolerance are structurally identical. A value-independent
("robust flat window") tolerance is just a 2-node spline with identical
values at both nodes — there is no separate flat-tolerance type; empirically
this is currently the only shape actually fit (see
`git/featureprediction`'s `AI.md`), value-dependence is supported but unused
so far. `rt_tol_sec`'s spline **values** are in seconds (converted to minutes
once at `Search` build time, `/60.0` — `scan_start_time` is in minutes); its
grid **x-axis** stays in whatever units the anchor points were fit on
(minutes, matching `scan_start_time` directly). `mobility_tol` needs no unit
conversion (1/K0 throughout).

The standalone `dump_peptides` binary (`crates/sage-cli/src/bin/dump_peptides.rs`)
digests a FASTA into a `peptide,proteins,monoisotopic,decoy` parquet without
needing spectra input or the full search config — only `--fasta` plus an
optional `database`-shaped JSON (enzyme/mods/mass bounds/decoy_tag) — for
`git/featureprediction` (a separate repo/pipeline stage) to generate
RT/IIM predictions from.

## Test fixture

`crates/sage-cloudpath/tests/data/pmsms_fixture/` is a small (~400K, 10-row)
real pmsms/precursors pair (both `.parquet` and `.mmappet` precursors
versions) used by `pmsms.rs`'s unit tests to assert both precursor formats
parse to identical spectra. `pmsms.mmappet`'s `mz` column (index 3, added
alongside the pre-existing `tof`/`intensity`/`score`) was back-filled by
applying the fixture's old `tof2mz.mmappet` lookup table once, offline, to its
`tof` column — same numeric values the old tof2mz-backed reader would have
produced, so this fixture's removal doesn't change what the unit tests assert.
**This repo's `.gitignore` has blanket `data/` and `*.txt` rules** that would
otherwise silently drop this fixture (the `.gitignore` has an explicit
negation carve-out for `crates/sage-cloudpath/tests/data/` — if you add
fixtures elsewhere under a directory named `data` or with
`.txt`/`.json`/`.csv`/`.tsv` files, check `git status`/`git check-ignore -v`
before assuming they're tracked).

## Testing

```
cargo test -p sage-cloudpath --lib                     # default features (no parquet)
cargo test -p sage-cloudpath --lib --features parquet  # matches the real `sage` binary build
```

`sage-cli` always builds `sage-cloudpath` with `features = ["parquet"]`
(see `crates/sage-cli/Cargo.toml`), so the actual `sage` binary always has
pmsms support; the `parquet`-feature-off path only matters if `sage-cloudpath`
is ever used as a library without it.
