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

RT and IIM are two fully independent, separately-optional filtering
dimensions — either, both, or neither may be configured. `--predicted-rt
<path.parquet>` (columns `sequence, rt`) loads a read-only
`HashMap<String, f32>`; `--predicted-iim <path.parquet>` (columns `sequence,
charge, iim`) loads a read-only `HashMap<(String, u8), f32>` — both once at
startup (`Runner::predicted_rt`/`predicted_iim`,
`crates/sage-cli/src/runner.rs`), shared across the parallel search same as
`database`. `rt_tol_sec` is required exactly when `--predicted-rt` is given
(and vice versa); same independent pairing for `mobility_tol`/
`--predicted-iim` — `Input::build()` rejects either half being set without
its partner, with no cross-coupling between the two pairs.
`Scorer::evict_rt_iim_mismatches` (`crates/sage/src/scoring.rs`) evicts
fragment-matched candidates whose predicted RT and/or IIM falls outside the
observed spectrum's window, right before `trim_hits` — each dimension
checked only if its map is configured at all; a candidate with no
`sequence`/`(sequence, charge)` entry in a configured map is left alone on
that dimension (permissive — avoids rejecting due to prediction-coverage
gaps). Design/history: `necromerge2`'s `plans/better_sage_filtering.md`
(original combined design) and `plans/rt_iim_independent_dimensions.md`
(the independent-dimension split, 2026-08-21 — replaced the single
`--predicted-properties`/`(rt, iim)`-tuple design outright, no back-compat
shim, since this fork has no external users).

Motivation for the split: a real F9477 comparison found RT and IIM filtering
don't contribute equally (Chronologer_RT beats SAGE's own internal RT model
by ~2.5x on robust residual error, while IM2Deep and SAGE's own IIM model
are much closer) — independent dimensions let that be measured directly at
the FDR level, and let RT-only filtering skip IM2Deep's Koina call entirely
(a real ~57-minute-per-run cost on a full human proteome, not just an
ablation nicety).

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

### Eviction lookup: dense peptide-index arrays, not string-keyed maps

`--predicted-rt`/`--predicted-iim` load `HashMap<String, f32>`/
`HashMap<(String, u8), f32>` from parquet (`Runner::predicted_rt`/
`predicted_iim`), but `Scorer` never reads those directly — `Runner::run`
resolves each into a dense, peptide-index-keyed `Vec<Option<f32>>` exactly
once, before the per-file search loop starts (`resolve_predicted_rt`/
`resolve_predicted_iim`, `crates/sage-cli/src/runner.rs`). This exists
because `evict_rt_iim_mismatches` originally called `Peptide::to_string()`
(itself doing a `unimod::lookup_reverse` per modified residue) twice per
fragment-matched candidate, on every spectrum — measured on a real F9477
run, combined RT+IIM eviction this way added +37% (152s) to `run_sage`'s
wall time versus no eviction at all. Resolving once, by peptide index
instead of by re-deriving each peptide's display string on every
(candidate, spectrum) pair, cuts that to a plain array read
(`self.predicted_rt.and_then(|by_idx| by_idx[peptide_idx])`) — real F9477
RT-only measurement: 473.4s (unoptimized) → 386.1s (optimized), ~18-20%
faster, within ~20s of the 406.6s no-eviction baseline.

IIM needs one extra step since it has a charge dimension `predicted_rt`
doesn't (Chronologer_RT has no charge input at all; IM2Deep requires one)
— `Scorer::iim_dense_slot` (`crates/sage/src/scoring.rs`) maps
`(peptide_idx, charge)` to a slot in a flat array of length
`n_peptides * (max_charge - min_charge + 1)`, shared between the build
side (`resolve_predicted_iim`) and the read side
(`evict_rt_iim_mismatches`) so the indexing arithmetic can't drift apart
between the two. This was **not** the first design tried — `predicted_iim`
was originally `HashMap<(usize, u8), f32>` (index-keyed, but still a hash
map). A standalone benchmark outside this crate (fixed-seed xorshift64*
PRNG for reproducibility, real F9477 scale — 18,902,646 entries, 20M
random-access queries, 80/20 hit/miss mix) measured:

| structure | ns/query | vs `std::HashMap` |
|---|---|---|
| `std::HashMap` (SipHash) | 239.1 | 1x |
| `FnvHashMap` (already a `sage` dependency, used elsewhere for the same reason) | 89.9-99.9 | ~2.4-2.7x |
| `FxHashMap` (rustc-hash) | 45.5-46.3 | ~5.2x |
| dense `Vec<Option<f32>>` | 10.4-10.8 | ~22-23x |

Dense wins outright — no hashing at all, and smaller (no per-entry key
storage or hashmap load-factor overhead: ~151MB vs the HashMap's
real-measured heavier footprint at this scale). A same-benchmark follow-up
comparing `Vec<Option<f32>>` against `Vec<f32>` with a negative sentinel
(both RT and IIM are always >= 0 physically) found **no query-latency
difference at all** (10.84ns both ways) — random access at this scale is
memory-*latency*-bound (one cache-line fetch per lookup, ~64 bytes,
regardless of the 4-vs-8-byte payload), not bandwidth-bound, so halving
the per-slot size doesn't remove any round-trips. `Vec<f32>`+sentinel
would still halve memory and build ~40% faster, but `Option<f32>` was kept
for the self-documenting safety (no sentinel-value convention to remember
or get wrong at a call site) since RAM headroom was never the actual
constraint on real hardware (~4-5GB peak for this whole feature on a
62GB/47GB-available machine).

`Runner::run` (which owns `self`) takes the raw string-keyed maps via
`.take()`, not `.as_ref()`, when resolving — the owned map moves into the
resolution closure and is dropped there, freeing it before the long
per-file search loop runs, rather than keeping both the raw and resolved
forms alive for the rest of the run. `prefilter_peptides` (the
`database.prefilter` path, off by default, no real job in this project
enables it) still uses `.as_ref()` since it re-resolves fresh per fasta
chunk and must keep borrowing across iterations — its own resolved arrays
are scoped to each chunk's short-lived mini database, never promoted to
`Runner` fields.

### External RT/IIM as LDA features, alongside the internal model (2026-08-24)

`--predicted-rt`/`--predicted-iim` above only ever fed a **hard** pre-hyperscore
eviction filter (`evict_rt_iim_mismatches`) — never SAGE's own `discriminant_score`
(the semi-supervised Fisher LDA in `crates/sage/src/ml/linear_discriminant.rs`,
used for FDR/q-values, not raw hyperscore). That LDA already had `sqrt(delta_rt_model)`/
`sqrt(delta_ims_model)` features, but those come from SAGE's own separate, much
weaker, in-run-only composition-regression model (`ml/retention_model.rs`/
`ml/mobility_model.rs`, fit on that run's own confident hits, called unconditionally
in `runner.rs` after search) — the external, per-run-calibrated predictions never
reached the LDA. Design/history: `necromerge2`'s `plans/lda_external_rt_iim_features.md`.

Closed by adding **two more** LDA features (`z2_rt_external`/`z2_ims_external`,
`((observed - predicted_external) / sigma) ** 2`) computed in `build_features`
from the same dense arrays `evict_rt_iim_mismatches` already reads, plus a new
`rt_sigma_sec`/`iim_sigma` config pair (robust MAD-based scale, from
`git/featureprediction`'s `tolerance.py::robust_sigma`, required alongside
`predicted_rt`+`rt_tol_sec` / `predicted_iim`+`mobility_tol` — three-way all-or-none
validation in `Input::build`). The internal model's own features are **kept**, not
replaced — this is additive, not a swap.

**Dynamic LDA column count, not a fixed-size slot with a `0.0` default.** SAGE's
LDA feature array was a `const FEATURES: usize = 20` fixed array
(`BASE_FEATURES`/`BASE_FEATURE_NAMES` now); a constant-valued extra column on any
run without external predictions configured would have zero within-class variance,
risking a singular covariance matrix in `LinearDiscriminantAnalysis::train` (`Gauss::solve`
returning `None`) — i.e. `discriminant_score` silently going uncomputed for the
*entire* run, not just a degraded feature. `score_psms` now takes `has_external_rt`/
`has_external_iim: bool` (from `Search.predicted_rt.is_some()`/`predicted_iim.is_some()`
in `runner.rs`) and builds a `Vec<f64>` per PSM instead of a fixed array, appending
the z² columns only when configured.

**z², not `sqrt(delta)` like the internal model's features.** The external
(Chronologer/IM2Deep) residual is already close to Gaussian after
`spectrum_q`-filtering (see `git/featureprediction`'s `confident_hits.py`), so z² is
already chi-square(1)-shaped — a reasonable LDA input as-is, unlike the internal
model's more skewed raw abs-residual (hence its `sqrt()` transform).

New `Feature` fields (`predicted_rt_external`, `delta_rt_z2_external`,
`predicted_ims_external`, `delta_ims_z2_external`) are output in both
`results.sage.tsv` and `results.sage.pin` (`runner.rs`) and the parquet schema
(`sage-cloudpath/src/parquet.rs`) — available for `sagepy-rescore`/mokapot to pick
up later, though nothing downstream consumes them yet (deliberately deferred, see
the plan doc).

**Explicitly deferred, not done here:** k-fold train/validation splitting of SAGE's
own LDA (mokapot-style) — `LinearDiscriminantAnalysis::train`/`score_psms` still
fit and score the same PSM set with no CV, same as before this change. Acceptable
for now given the LDA's low parameter count (20-22 coefficients) relative to typical
PSM counts; revisit separately if it becomes a concern.

## Unimod modification support (`crates/sage/src/unimod.rs`)

Found (2026-08-21) empirically against the live Koina server: both
`Chronologer_RT` and `IM2Deep` silently ignore this fork's own `[+mass]`
mass-delta bracket notation for modifications — `PEPTC[+57.0216]DEK` and
`PEPTCDEK` come back with byte-identical predictions, only genuine
`[UNIMOD:<id>]` notation actually changes anything. Traced to Koina's own
`IM2Deep_Preprocess_AC`/`Deeplc_Preprocess_onehot` server-side code
(`wilhelm-lab/koina` GitHub): internal mods are silently regex-skipped
(no error), N-terminal mods used to crash outright (`KeyError: '+42'` —
this fork's own N-term notation, `[+mass]-SEQUENCE` with a hyphen, *does*
get split correctly as a terminal-mod segment, then the raw mass string
gets looked up as if it were itself a literal Unimod ID).

`static_mods`/`variable_mods` (`database` config) now accept a modification
value as either a plain float (unchanged) or a `"UNIMOD:<id>"` string,
resolved to that id's exact monoisotopic mass at config-load time
(`crates/sage-cli/src/input.rs`'s `resolve_unimod_refs`/
`resolve_unimod_refs_in_database` — walks the raw config JSON *before*
deserializing into `Builder`, so `Builder.static_mods: HashMap<String, f32>`
itself never changes shape; the rest of the engine only ever sees plain
floats, exactly as before). A plain float that coincides (within 5 mDa —
`unimod::COINCIDENCE_TOLERANCE_DA`) with a real Unimod entry is a config
**error**, not silently accepted — e.g. `57.0216` (what real job configs
actually write for Carbamidomethyl) errors naming `UNIMOD:4`, since letting
it through would leave that modification permanently unable to round-trip
as `[UNIMOD:4]` on output despite plainly being that modification.

Resolution is against an *active* table — `crates/sage/data/unimod.csv`
(`id,name,mono_mass`, regenerate via `scripts/build_unimod_table.py`,
same source `https://www.unimod.org/obo/unimod.obo` Koina's own
preprocessing downloads) embedded into the binary via `include_str!` and
used by default, or `--unimod-db-path <path>` (both the `sage` binary and
`dump_peptides`) to override it entirely with an external CSV of the same
shape. Either way, parsed once, lazily, into a process-global
`OnceLock` — correct for the real one-shot-CLI-process use case; test
code must take care not to have more than one test per binary actually
resolve a `UNIMOD:<id>` reference or override the active table (same
constraint, and same reasoning, as the existing capturing-logger tests in
`crates/sage-cli/src/input.rs`).

`Peptide::Display` (`crates/sage/src/peptide.rs`) round-trips: a
modification whose mass is bit-identical to one resolved via `UNIMOD:<id>`
*this run* prints `[UNIMOD:<id>]` instead of `[+mass]` — deliberately
provenance-scoped to just the ids this run's config actually referenced
(`unimod::set_reverse_table`/`lookup_reverse`), not a full ~1560-entry
reverse lookup, so this never fires on a coincidental match the config
didn't explicitly ask for. `dump_peptides`' own output, `results.sage.tsv`'s
`peptide` column, and this fork's own `--predicted-properties` `HashMap`
lookup key (`scoring.rs`) all stay consistent with each other automatically
as a result — no call-site changes needed anywhere `Peptide::to_string()`
was already used. When no `UNIMOD:<id>` reference is used at all, output is
byte-for-byte identical to before this existed.

Verified against the real, live Koina server that exposed the bug (not
just unit tests): real `dump_peptides` output for
`[UNIMOD:1]-MPEPTC[UNIMOD:4]DEK` (N-terminal + internal mod, the exact
shape that used to crash) round-trips through both `Chronologer_RT` and
`IM2Deep` without error, and the modified/unmodified predictions now
genuinely differ (e.g. IM2Deep CCS 338.67 → 346.43 Å² for the internal
mod alone).

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
