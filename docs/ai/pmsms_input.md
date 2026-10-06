# pmsms/precursors input

Core reader: `crates/sage-cloudpath/src/pmsms.rs`. Inputs, given as explicit
paths (never assumed to live together in one fixed-filename directory):

- `pmsms.mmappet/` — schema.txt + one `<idx>.bin` file per column (see
  `git/mmappet`'s Python writer), columns looked up by name. Needs
  `tof`(uint32) and `intensity`(uint32). An `mz` column, if present, is
  ignored with a warning.
- `tof2mz` table — an mmappet with one float32 or float64 column `mz` indexed by
  tof, read as f64; fragment m/z is `table[tof]`, rounded to f32. A tof outside the
  table is an error (`PmsmsError::TofOutOfRange`).
- precursors, either `.parquet` (via the `parquet` crate) or `.mmappet/`
  (schema.txt + one `<idx>.bin` file per column, columns looked up by name —
  see `git/mmappet`'s Python writer for the on-disk format). Needs 7 columns:
  `precursor_idx`(u64), `mz`(f64), `rt`(f64), `inv_ion_mobility`(f64),
  `charges`(i64), `fragment_spectrum_start`(u64), `fragment_event_cnt`(u64).
- optionally a per-precursor `fragment_shift_ppm`(f64) column (below).

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

Any pmsms read error (a missing or wrong-dtype column, a tof outside the table, an
unreadable precursors table) panics in the runner, so SAGE exits non-zero. Until
2026-10-05 it was logged and turned into an empty spectrum list: SAGE searched 0
spectra, wrote valid-looking empty outputs and exited 0, and the failure surfaced only
downstream.

## Fragment m/z recalibration: a recalibrated table plus `fragment_shift_ppm` (2026-10)

The fragment correction fitted in necromerge2 (searchops' `recalibrate_pmsms_mz`) is
`bias + f_mz(mz) + f_rt(precursor rt)` in ppm. SAGE gets it in two pieces, neither
of which needs any model code here:

- `f_mz` depends only on m/z, which depends only on tof, so it is folded into a
  float64 tof2mz table, `table[t] / (1 + f_mz(table[t])·1e-6)`, passed as `--tof2mz`.
- `bias + f_rt(rt)` is constant per precursor and arrives as the precursors column
  `fragment_shift_ppm`. When the column is present, each precursor's fragment m/z
  are `(table[tof] / (1 + fragment_shift_ppm·1e-6)) as f32`, computed in f64
  (`look_up_mz`); when absent, `table[tof] as f32`. No flag: like
  `ppm_tol_lo`/`ppm_tol_hi`, the column's presence is the switch, so nothing checks
  that a shifted precursors table is paired with the recalibrated table — the
  pipeline wires both together.

The table is float64 because a float32 one would round each m/z twice: measured on
F9477, that moves 25% of peaks by one float32 ulp against a single rounding, while
float64 moves essentially none. With the raw float32 table and no column the
arithmetic (`f32 -> f64 / 1.0 -> f32`) returns the table value exactly.

The mzML writer (`pmsms2mzml --tof2mz`), the MGF writer (`msms2mgf --tof2mz`) and
`sage-pmsms-mapper --tof2mz` perform the same float64 operations from `tof`, so exports
hold exactly the m/z SAGE searched with no materialized `mz` column (necromerge2
`plans/exports_from_tof2mz_table.md`).

Why the per-precursor part is computed upstream from raw RT rather than by SAGE: in
necromerge2's mode 3 the RT correction overwrites the precursors' `rt` column before
the final search, while `f_rt` was fitted on raw RT. The column also leaves the
per-precursor model free to change (e.g. a 2-D `f(rt, 1/K0)`, or a free intercept)
without touching SAGE.

History: the first version (git/sage `6e06d5d`) took the `f_mz` curve itself,
`--mz-recalibration <dir.mzcalib>`, and evaluated it per peak with an f64 copy of
timstofu's `EvenlySpacedLinearSpline.njit_evaluator`, pinned bit for bit against
numba-generated values, plus `.mzcalib` dimension checks. Folding `f_mz` into the
table removed the flag, the spline, the `.mzcalib` reader and its error cases.

## Matched peaks' pmsms rows: `fragment_pmsms_row`, `closest_fragment_pmsms_row` (2026-10)

The matched-fragments export (TSV and parquet) carries, for the matched (most intense)
peak and the closest peak, the row of that peak in the input pmsms, so any pmsms column
is one lookup away. `-1` when the input is not a pmsms.

Processing discards a peak's input position (deisotoping, top-N selection, the mass
sort), so the row is recovered after scoring instead of carried through it:
`score_candidate` keeps each matched/closest peak's processed mass and charge scale in
`Fragments`, and the runner's `resolve_pmsms_rows` asks `pmsms::PmsmsRowResolver`, which
re-opens the pmsms `tof` column, the tof2mz table and the precursors, recomputes the
masses of the precursor's slice with the reader's own float32 arithmetic
(`((table[tof] / divisor) as f32 - PROTON) * charge`) and binary-searches it (the slice
is tof-sorted, hence mass-sorted; a linear scan backs it up). Equal float32 masses
resolve to the lowest row; a deisotoped peak resolves to its monoisotopic root, whose
m/z deisotoping keeps.

F9477: all 2,245,178 matched fragments resolved, none to `-1`, the m/z at every row
reproducing `fragment_mz_experimental` bit for bit; 354 ms for the whole export. It
replaced searchops' `sage-pmsms-mapper` (nearest m/z within 0.001 Da), which found the
right rows per precursor but, sorting each PSM's m/z without its fragment indices,
attributed 89% of them to the wrong SAGE fragment.

## Borrowed fragment input and processed storage

The explicit-path CLI input calls `parse_processed`: each precursor's m/z
(looked up as `table[tof]`, divided by the precursor's shift, into a
per-spectrum vector) and mapped `intensity: &[u32]` range feed `SpectrumProcessor` directly,
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
