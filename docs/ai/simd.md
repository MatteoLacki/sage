# SIMD fragment-page search on F9477

Branch: `optimizations/simd`, based on `42510e4` (`Add unlimited_max_peaks and assume_sorted_peaks config flags`).

## Correctness update, 2026-09-15

Both original benchmark binaries contain the pre-existing per-kind fragment-index
regression introduced by `5753d07`. Their output equivalence establishes that the
index optimization preserved the affected behavior, not that MS2 similarity
features were correct. Fix `548ffd5` restores per-kind numbering. See
[`interai/2026-09-15_ms2_similarity_regression.md`](../../interai/2026-09-15_ms2_similarity_regression.md).
The original timings below remain historical measurements; spectrum-layout and
mapped-input follow-ups must use a corrected `548ffd5` control and corrected
optimized builds. The isolated one-ULP posterior-error variation is separate
from this systematic fragment-annotation and intensity-lookup regression.

## Change

The fragment index stores separate `f32` masses and `PeptideIx` columns, preserving the original fragment-page and within-page peptide order. Its resident payload remains eight bytes per fragment. Index construction still sorts temporary `Theoretical` records before splitting them into the two columns, so construction temporarily holds both representations.

`page_search_batch` resolves an exact, half-open precursor-compatible peptide-ID range, searches each page's ID column, then scans only the corresponding mass slice. On x86-64 with runtime-detected AVX2, slices of at least eight elements use eight-lane ordered mass comparisons and a match mask. Matching lanes are emitted in ascending order; peptide counter updates remain scalar. Short slices, tails, and CPUs without AVX2 use the scalar implementation.

Duplicate fragments and overlapping/repeated windows retain their multiplicity. No scoring arithmetic, tolerances, search configuration, or floating-point compiler flags change. The public Rust storage API changes from `fragments: Vec<Theoretical>` to `fragment_mzs` / `fragment_peptide_indices`; `fragments()` reconstructs records lazily. Search callbacks/iterators now yield these small records by value. Workspace callers are updated.

This experiment measures the combined layout, boundary handling, and SIMD change. It does not isolate SIMD's contribution from the other two changes.

## Dataset and method

- Full F9477: 716,614 precursor rows, human FASTA, 6,300,882 peptides and 198,932,246 indexed fragments.
- Source pipeline node: `nodes/run_sage/19b746debf1e961d99207116c6132e96c259413eca5d6d1505ce11c90673273d`. Its dependency graph traces to `source_bruker_d/6fd11fbc1a5797a93fd57b1e4db6f6dd2afcb8543c194ae1cceac7e7861dba7f`, configured with `data/F9477.d`.
- Existing production configuration: bucket size 32,768; max 800 peaks; fragment charge 1; isotope errors 0–3; deisotoping; RT filtering and predicted fragment-intensity features; PIN and fragment annotations enabled.
- AMD Ryzen 7 4800H, 8 cores / 16 threads. `RAYON_NUM_THREADS=16`. Normal release profile, fat LTO, one codegen unit, no additional Rust flags.
- Both binaries use identical spectra, precursor table, FASTA, configuration and prediction values. Outputs go only to `target/simd-benchmark`. Telemetry disabled.
- Existing predictions had the old sequence-keyed schema. `prepare.py` converts local copies to the current positional schema, checks complete RT coverage, matches fragment entries by sequence/charge, preserves missing-entry sentinels, and embeds the peptide-mass fingerprint. Current Sage verifies that fingerprint against its own digest. Production nodes remain untouched. The initial rejected run (`baseline-1`) is excluded from timings.
- Search time comes from Sage's search-phase log. Wall time and peak RSS come from `/usr/bin/time -v`; the wrapper also records elapsed monotonic time. Runs are sequential; compilation, tests, and large analysis jobs do not overlap measured runs. Page caches are not dropped. First full pair and repeat pair are separated in time; compare paired runs as well as their aggregate.

## Measurements

| Run | Search (s) | Wall (s) | Peak RSS (GiB) | User CPU (s) | System CPU (s) |
|---|---:|---:|---:|---:|---:|
| baseline-2 | 315.719 | 371.754 | 28.861 | 4916.49 | 43.48 |
| simd-1 | 290.365 | 339.756 | 28.860 | 4496.82 | 36.09 |
| simd-2 | 284.234 | 343.212 | 28.903 | 4475.64 | 53.73 |
| baseline-3 | 336.457 | 388.010 | 28.839 | 5004.27 | 39.68 |

Mean search time: **326.088 → 287.300 s (11.90% reduction; 1.135× throughput)**. Mean wall time: **379.882 → 341.484 s (10.11% reduction)**. Peak RSS is effectively unchanged, approximately 28.9 GiB.

Paired search-time reductions are 8.03% and 15.52%. Only two runs per version were measured on a shared workstation, without CPU-frequency pinning; this is an observed improvement with visible run-to-run variation, not a precise universal speedup. User CPU time also falls, from a mean 4,960.38 s to 4,486.23 s (9.56%).

Compiler: rustc 1.97.1 (`8bab26f4f`, LLVM 22.1.6). Frozen binary SHA-256:

- Baseline: `d8e5adfaaf124de357e77a8a7846497692f4af612ba1e93efa1b8666112b45e5`
- Optimized: `4dea2823b7d968ce0d3ca149f9fb8af12c766aaebb9ee585fcaf57b76d564076`

## Validation

`cargo test --workspace --all-features --offline`: 200 tests passed. New tests compare page-search matches against a brute-force oracle across page sizes, exact/isobaric precursor boundaries, empty searches, repeated fragments and repeated windows. AVX2 and scalar scans are compared across unaligned starting positions, lengths covering every tail, exact boundaries, infinities, NaNs and reversed windows.

The release binary contains `vmovups`, `vcmpge_oqps`, `vcmple_oqps`, `vandps` and `vmovmskps` operating on YMM registers in `scan_masses_avx2`.

Both optimized runs have identical 363,156 PSM identities and all 59 non-ID TSV columns byte-identical to the first baseline, after matching/sorting by spectrum, peptide, label, charge, rank and isotope error. Scheduling-dependent `psm_id` differs.

The repeat baseline has one additional difference, also present when comparing the two **baseline** runs: `posterior_error` for `precursor_idx=616655`, peptide `QEYDESGPSIVHR`, changes from `-224.89215` to `-224.89214` (one float32 ULP). Both optimized runs retain the first baseline's value. All other non-ID columns match. This is existing downstream run-to-run variation, not a difference introduced by the patch. KDE density estimation (`ml/kde.rs`) uses a parallel floating-point fold and sum; changing reduction order is a likely mechanism, though that particular reduction was not instrumented. The strict comparison script intentionally flags it; its failure report is retained as `compare-baseline-3-simd-2.json`.

## Reproduction of the historical index experiment

Run from this repository root. Python preparation/comparison use the pipeline's existing `../../venvs/common/bin/python`; the run wrapper needs only Python's standard library. The scripts are specific to this F9477 fixture and write under `target/simd-benchmark`.

1. Preserve a release baseline binary built from `42510e4` as `target/simd-benchmark/sage-baseline`.
2. Run `../../venvs/common/bin/python docs/ai/simd_benchmark/prepare.py` to recover the saved invocation and prepare prediction copies.
3. Build commit `01600dd` with `cargo build --release --bin sage --offline`; preserve the result as `target/simd-benchmark/sage-simd`. Later branch revisions include additional changes and do not reproduce this historical binary.
4. Run `python3 docs/ai/simd_benchmark/run.py LABEL BINARY` for each timing. Labels must be new output directories. The original run sequence was baseline-2, simd-1, simd-2, baseline-3.
5. Run `../../venvs/common/bin/python docs/ai/simd_benchmark/compare.py BASELINE_LABEL SIMD_LABEL` for exact TSV comparisons.

`target/simd-benchmark/manifest.json` records actual input paths and conversion provenance. Each run directory contains `command.json`, `run.log`, `time.txt`, `summary.json`, and complete Sage outputs. The benchmark scripts and this report are tracked; large binaries, converted inputs and outputs are ignored build artifacts.

## SIMD utilization diagnostic

A separate instrumented build searched a uniform 10,000-of-716,614 precursor sample (`pandas.sample(random_state=9477)`), retaining the full database and production search configuration. Instrumentation counts precursor-compatible mass-slice lengths per fragment window, using per-worker counters. It is absent from the measured full-run binaries and restored out of the source after the diagnostic build.

- 37,269,240 window/range scans.
- 22,366,979 empty ranges (60.01%).
- 160,676 ranges contain at least eight masses (0.431%).
- 29,891,644 mass comparisons; 2,199,376 (7.36%) fall in complete eight-lane blocks.
- Mean precursor-compatible range length: 0.802 fragments.

Most narrow-window F9477 scans cannot benefit from eight-lane filtering. These counts suggest layout and exact-boundary handling are important contributors to the measured improvement; a scalar-SoA ablation would be needed to assign the speedup quantitatively. SIMD remains useful for longer ranges, whose prevalence depends on precursor tolerance and database distribution. The diagnostic's 4.129 s search time is not included in full-run comparisons.

To reproduce the diagnostic after the full measurements: run `sample.py` with the pipeline Python, `build_diagnostic.py` with Python 3, then `run_diagnostic.py` with the pipeline Python, all from `docs/ai/simd_benchmark/`. `build_diagnostic.py` saves/restores both affected Rust files in a `finally` block and restores the frozen optimized release binary. Avoid concurrent edits to those two files while it runs. Raw counts are retained in `target/simd-benchmark/scan-statistics.json` and the diagnostic log.

## Processed spectrum SoA and borrowed input, 2026-09-15

The corrected control is commit `548ffd5`: SIMD fragment index plus fixed
per-kind fragment numbering, with the original processed AoS/raw-copy path.
The SoA-only snapshot replaces non-mobility `Vec<Peak>` storage with
`PeakColumns`, searches the mass column directly, and passes the existing
intensity slice to top-six scoring rather than allocating a per-candidate copy.
Chimera removal compacts mass, intensity and observed charge columns together.

The final snapshot additionally processes borrowed mapped m/z and u32 intensity
ranges for explicit `--pmsms`/`--precursors` input. Parallel collection preserves
precursor order; returned spectra own their selected SoA columns. Conversion to
f32 precedes intensity comparisons and accumulation, preserving rounding ties.
Final columns reserve the exact retained count through an exact-size iterator;
temporary selection/deisotoping storage still scales with raw peak count.
Legacy owned-reader APIs and positional `.pmsms` input still copy raw columns.
See [input details](pmsms_input.md).

Both optimized snapshots retain fix `548ffd5`. The new scoring regression test
uses observed spectra matching nonuniform predictions, checks both b/y traversal
orders, exact fragment ordinals, and unit cosine/entropy/Pearson similarity.
Temporarily restoring global numbering makes this test fail on zero/negative
y-ordinals; the fixed source was restored before builds and measurements.
Borrowed-input tests additionally cover u32-to-f32 rounding ties, real isotope
envelopes, empty/singleton spectra, top-N limits, and both precursor formats,
including metadata and use of returned spectra after mapping closure.

### Measurements

Same full F9477 fixture and 16-thread release settings as above. Each row is
one complete run; the last control repeats the frozen corrected-control binary.
The first control occurred several hours before the other runs, so the repeat
helps expose workstation/cache drift. No compilation, tests, or large output
comparisons overlap measured runs.

| Run | Input + preprocessing (s) | Search (s) | Wall (s) | Peak RSS (GiB) |
|---|---:|---:|---:|---:|
| fixed-control-1 | 31.945 | 286.132 | 344.456 | 28.868 |
| fixed-soa-2 | 25.213 | 288.222 | 339.375 | 28.856 |
| borrowed-2 | 14.270 | 291.286 | 330.353 | 20.232 |
| fixed-control-2 | 27.559 | 352.777 | 407.408 | 28.858 |

The robust benefit is memory and input/preprocessing: peak RSS falls from
approximately 28.87 to 20.23 GiB
(29.9% lower). Input plus preprocessing takes
14.270 s with borrowing versus 27.559 s
in the repeat control (48.2% lower), and 25.213 s
in the SoA-only run. SoA alone does not establish a search-phase gain; one run
per optimized variant and an unpinned shared workstation limit timing conclusions.
These follow-ups measure changes on top of the SIMD index, not a remeasurement
of the historical index speedup. No scalar-SoA index ablation was run.

### Validation and artifacts

`cargo test --workspace --all-features --offline`: **205 tests passed** for the
final source. SoA-only had 203 passing tests. Both optimized variants have the
same 363,156 PSM identities as corrected control and all 59 non-ID TSV columns
byte-identical, including every MS2 similarity feature and posterior error in
these comparisons. The final borrowed-input build additionally matches all
2,282,985 fragment annotation rows exactly after PSM-ID mapping and row sorting;
every ordinal lies between 1 and peptide length minus one.

Frozen binaries and SHA-256:

- `sage-fixed-control`: `2f495ec580f6f843b39ea2b6e9ccc43c03f16c7fcbcbb970af013976934d930f`
- `sage-fixed-soa`: `6be11d780f82c8e953eacf61dad96670008834f0cafb70370775f4db30fe5ab5`
- `sage-borrowed`: `9a9ac63ff6ca38fec97a337e3d22a1e11ffed074f254fa00e489b8b438366e56`

Build snapshots are `fixed-soa.patch` and `borrowed.patch` against `548ffd5`,
under `target/simd-benchmark`. The corrected control was built from an isolated
`git archive 548ffd5`. When sharing Cargo's target directory between source
snapshots, force a rebuild of the three local packages before freezing the next
executable; a cached build can otherwise leave the previous executable at the
shared `target/release/sage` path. Verify distinct binary hashes before comparison.
Separate target directories also avoid that collision.

Excluded attempts: `fixed-soa-1` was stopped after hash verification identified
the control executable at the shared output path; `borrowed-1` stopped before
launch because its build was not yet ready. Neither contributes measurements.
The run wrapper now records each binary's SHA-256 and input/preprocessing time.

`target/simd-benchmark/spectra-timings.json` preserves the full table and CPU
measurements. Regression-test pass/mutation logs, workspace-test logs, build
logs, full outputs, and comparison JSON files are retained alongside it.
Reproduce PSM comparison using `compare.py LEFT RIGHT`; use
`compare_fragments.py LEFT RIGHT` for exact annotations and ordinal bounds,
both with the pipeline Python from this repository root. Run comparisons after
measurements, since sorting full annotation tables consumes CPU and memory.
