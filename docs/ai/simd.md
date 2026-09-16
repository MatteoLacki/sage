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

## Interleaved bound resolution and `bucket_size`, 2026-09-15

Starting point `6c20f3f`. The control binary rebuilt from that commit hashes
`9a9ac63f…`, identical to the frozen `sage-borrowed` above, so the two
measurement rounds share a control.

### What the previous round left unmeasured

The SIMD utilization diagnostic counted scans and lane occupancy but never
attributed elapsed time, so "most narrow-window scans cannot benefit from
eight-lane filtering" was never turned into a bound on what the kernel could
be worth. Combining that diagnostic's own two outputs does turn it into one:
it spent 4.129 s × 16 threads = 66 CPU-s on 29,891,644 mass comparisons, i.e.
2.2 µs per comparison. A float compare is ~1 ns, so the comparisons were never
where the time went, and no change to the comparison kernel could have returned
more than a few percent.

A second instrumented build (`instrument_latency.py`) confirmed this directly by
bracketing the two phases of `page_search_batch` with `lfence`-fenced `rdtsc`
reads, on the same 10,000-precursor sample:

| Phase | Cycles | Share |
|---|---:|---:|
| Peptide-ID bound resolution | 170,703,174,353 | **94.1%** |
| Fragment-mass scan (incl. AVX2 kernel) | 10,632,892,021 | 5.9% |

28,366,868 page groups, 37,269,240 windows — 2,837 groups per spectrum, each
paying two `partition_point` calls over the DRAM-resident peptide-ID column.
That is 6,017 cycles per group, ~200 cycles per binary-search step: full memory
latency, with nothing overlapping it. The instrumented run took 4.374 s against
4.129 s clean, so the 5.9% instrumentation overhead does not distort the split.

The fragment-mass scan is the *only* phase the existing AVX2 kernel touches.
Its 5.9% share caps that kernel's whole addressable surface at a 1.06× speedup.

### Change

`page_search_batch` now splits into three phases: collect distinct pages, resolve
every page's `[inner_left, inner_right)` in `resolve_page_peptide_bounds`, then
scan. The resolver steps 16 pages' searches in lockstep, two independent searches
each (lower and upper bound), so 32 cache misses are outstanding at once instead
of one. `inner_right` is resolved against the whole page rather than
`ids[inner_left..]` — the column is sorted ascending and `pre_idx_lo <=
pre_idx_hi`, so the index is the same, and decoupling them doubles the lanes.
The inner step is branchless with a fixed trip count; a converged lane gets
`half == 0` and re-probes its own `base` without moving.

A binary search is a dependent load chain — the next address comes from the value
just loaded — so a single search leaves the core's ~22 outstanding-miss slots idle
while it stalls. Interleaving is what fills them. This is a memory-parallelism
change, not a vectorization one: AVX2 gather would express the same probes in one
instruction but Zen2 decomposes gathers into individual loads, so it would not
beat the scalar interleaving.

`bucket_size` turned out to matter more than the code. It was swept as a
config-only experiment, no rebuild:

| `bucket_size` | Search (10k sample, optimized build) |
|---:|---:|
| 8,192 | 2653 ms |
| 32,768 (production) | 1771 ms |
| 65,536 | 1573 ms |
| 131,072 | 1306 ms |
| 524,288 | 755 ms |
| 2,097,152 | 451 ms |
| 4,194,304 | 411 / 390 ms |

Larger pages are faster because page-group count, not page length, sets the probe
budget. Fragment windows are narrow (±13 ppm ≈ 0.026 Da) and mostly land on one
page whatever the page size, so 128× larger pages cut groups per spectrum from
2,837 to under 48 — 23.8× fewer binary searches — while each search grows by only
seven steps. The precursor-filtered mass slice grows in proportion, but that slice
is scanned sequentially and prefetchably, and is exactly what the AVX2 kernel
accelerates. Random dependent probes are traded for streaming scans.

### Measurements

10,000-precursor sample, isolating the two contributions:

| Build | `bucket_size` | Search | vs production |
|---|---:|---:|---:|
| control `9a9ac63f` | 32,768 | 4071 / 3966 ms | 1.00× |
| optimized `cd2c3da5` | 32,768 | 1771 ms | 2.27× |
| control `9a9ac63f` | 4,194,304 | 534 ms | 7.53× |
| optimized `cd2c3da5` | 4,194,304 | 411 / 390 ms | **10.0×** |

`bucket_size` alone accounts for 7.53×; the interleaved resolver adds 1.37× on
top of it, and is worth 2.27× on its own at the production page size.

Full F9477 (716,614 precursors), paired runs, same session and load, 16 threads:

| Run | Search (s) | Wall (s) | User CPU (s) | Peak RSS (GiB) |
|---|---:|---:|---:|---:|
| control `9a9ac63f`, bucket 32,768 | 279.083 | 317.72 | 4446.22 | 20.23 |
| optimized `cd2c3da5`, bucket 4,194,304 | **29.791** | **68.22** | **731.30** | 20.23 |

**Search 9.37×, wall 4.66×, user CPU 6.08×.** Peak RSS is unchanged — this buys
time with neither memory nor accuracy.

### Where the time goes now

Re-running the phase instrumentation on the optimized build shows the balance has
inverted:

| | bucket 32,768 | bucket 4,194,304 |
|---|---:|---:|
| Page groups | 28,366,868 | 1,192,620 |
| Bound-resolution cycles | 57,564,771,770 (89.7%) | 2,913,117,889 (41.4%) |
| Mass-scan cycles | 6,622,371,531 (10.3%) | 4,129,274,968 (58.6%) |

Two consequences. The resolver cut bound-resolution cycles 2.97× at the production
page size (170.7e9 → 57.6e9), confirming the mechanism rather than just the wall
clock. And at the large page size the mass scan is now the *majority* of
`page_search_batch` — so the AVX2 kernel, capped at 1.06× when this round started,
is now operating on the dominant phase. Any further work on that kernel should be
measured in this regime, not the old one.

`page_search_batch` itself is now only ~39% of search CPU (7.04e9 cycles vs
435 ms × 16 threads), against ~87% before. The remaining ~61% is scoring, which
was never profiled and is the next thing to attribute.

### Caveats

`bucket_size = 4,194,304` leaves 48 pages for this database. That is near-degenerate
— the mass bucketing is nearly vestigial — and the optimum depends on database size,
fragment tolerance and peak count, none of which were varied. Returns are already
flat from 2,097,152 (451 ms) to 4,194,304 (411 ms), so an intermediate value in the
1M–4M range captures nearly all of the gain with a less extreme structure. A value
chosen from index size rather than hardcoded would be better than any constant here.

Only one or two runs per configuration, on an unpinned shared workstation. The
9.37× full-run figure comes from a single paired comparison; the effect is far
larger than the run-to-run variation seen in earlier rounds (which reached 15%),
but it is not a precisely bounded number.

### Validation

`cargo test --workspace --all-features --offline`: **206 tests passed** (205 before,
plus `batched_page_bounds_match_partition_point`, which checks the interleaved
resolver against the per-page `partition_point` pair it replaced across bucket
sizes 1–17, one- to nine-page indices, partial final pages, duplicate IDs, and
every lower/upper key pair including empty and full-page ranges).

Full F9477 output equivalence against the paired control:

- 363,156 PSM identities identical; all 59 non-ID TSV columns byte-identical.
- All 2,282,985 fragment annotation rows exact after PSM-ID mapping and row
  sorting; every ordinal within 1..peptide length − 1.
- Only `psm_id` differs, as in every previous round (scheduling-dependent).

Additionally checked at the byte level, rather than through the comparison
scripts' per-cell string equality: dropping the scheduling-dependent ID column
and sorting the remaining lines makes `results.sage.tsv`,
`matched_fragments.sage.tsv` **and `results.sage.pin`** hash-identical between
the two runs. `results.sage.pin` is covered by neither `compare.py` nor
`compare_fragments.py` — it is a separate output with its own feature columns,
and its ID column is `SpecId`, not `psm_id`. Any future round should check it
explicitly; the comparison scripts alone do not.

Output equivalence was additionally confirmed on the 10k sample at bucket sizes
524,288, 2,097,152 and 4,194,304, establishing that `bucket_size` is purely an
index-organization parameter over a 128× range.

Binary SHA-256: control `9a9ac63ff6ca38fec97a337e3d22a1e11ffed074f254fa00e489b8b438366e56`,
optimized `cd2c3da59ed9f3296856a2e4c05a02049e25ae45473ed9a7340872f0a8702cb2`.

### Reproduction

`perf` was unavailable (`perf_event_paranoid = 4`, needs root), so all attribution
is in-source. `instrument_latency.py` / `build_latency_diagnostic.py` patch, build
and restore as `instrument.py` does; `run_sample.py LABEL BINARY [CONFIG]` and
`run.py LABEL BINARY [CONFIG]` run the 10k sample and the full fixture, both taking
an optional config path to override `bucket_size`. Config variants are written to
`target/simd-benchmark/config_bucket<N>.json`. Freeze each binary and check its
SHA-256 before comparing — the shared `target/release/sage` collision documented
above still applies.

## Cross-dataset replication on F9468, 2026-09-15

F9468 is an independent dataset, **1,708,243 precursors (2.4x F9477's 716,614)**,
with near-identical recalibrated tolerances: precursor ±6.08/6.43 ppm (F9477:
±6.57/6.76), fragment ±13.24/12.34 ppm (F9477: ±13.33/12.32), same `max_peaks:
800` and `max_fragment_charge: 1`. So it varies dataset *size* by 2.4x while
holding the search regime fixed.

Four frozen binaries, one per generation of this work, each at both page sizes.
Full dataset, 16 threads, sequential runs:

| Commit | Change | bucket 32,768 | bucket 4,194,304 | bucket gain |
|---|---|---:|---:|---:|
| `440c6f8` | before index sharing | 738.6 s | 269.3 s | 2.74x |
| `42510e4` | + `page_search_batch` sharing, pre-SIMD | 748.9 s | 180.6 s | 4.15x |
| `6c20f3f` | + columnar pages / AVX2 | 649.3 s | 87.4 s | 7.43x |
| `950641f` | + interleaved bound resolution | **293.9 s** | **65.4 s** | 4.49x |

Production-equivalent control (`6c20f3f` at 32,768) to best (`950641f` at
4,194,304): **9.93x**. Oldest measured commit to best: **11.29x**.

Every effect reproduces within a few percent of F9477:

| Effect | F9477 | F9468 |
|---|---:|---:|
| `bucket_size` alone, on `6c20f3f` | 7.53x | 7.43x |
| `bucket_size` alone, pre-SIMD `42510e4` | 4.24x | 4.15x |
| Interleaved bounds alone, at 32,768 | 2.27x | 2.21x |
| Overall control -> best | 9.37x | 9.93x |

The qualitative findings reproduce too: index sharing is **negative** at 32,768
(738.6 -> 748.9 s) and positive only once pages are large (1.49x at 4M); the
columnar/AVX2 change is worth 1.15x at 32,768 but 2.07x at 4M. Both are starved
by page geometry, in the same direction and for the same reason as on F9477.

### Validation

`f9468-mlp-control-b32768` against `f9468-mlp16-b4194304`: **640,621 PSM
identities** identical, all 59 non-ID TSV columns byte-identical, and all
**3,730,457** fragment annotation rows exact after PSM-ID mapping and row
sorting, every ordinal within 1..peptide length - 1. Only `psm_id` differs.

### Inputs and how they were produced

`jobs/f9468_benchmark.toml` in the parent `necromerge2` repo is a byte-for-byte
copy of `jobs/f9477_best.toml` with `tdf_path` swapped to `data/F9468.d`, so
recalibration settings are identical and the comparison is not confounded by
config drift. `./nf jobs/f9468_benchmark.toml -call` produced the data-dependent
nodes in ~15 min; the ~90 min of RT/IIM/fragment-intensity prediction did **not**
rerun, being peptide-derived and already cached for the same FASTA and database
config.

Two things to know when repeating this:

- The pipeline's own final `run_sage` node **fails**, unrelated to any change
  here: the cached `export_fragment_intensity_for_sage` parquet predates the
  positional-schema migration and lacks `dumped_peptides_sha256`, which current
  Sage rejects. The benchmark sidesteps it by reusing the already-prepared
  `target/simd-benchmark/fragment_index.parquet`, which carries the fingerprint.
  That substitution is sound because both peptide dumps
  (`dump_peptides/e3c5b424…` and `…/1fdcfee7…`) are byte-identical
  (`e90910feecb0e698e6804412b28d1066030d9d9c729e72e8f9edb0f124482b89`) --
  verify that before reusing it again. The node needs regenerating for the
  pipeline itself to complete.
- `source_dump_peptides_binary` is a symlink to the live
  `git/sage/target/release/dump_peptides`, so rebuilding that binary makes the
  node stale and gives `dump_peptides` a new hash. `ln -s` there is not
  idempotent: a pre-existing symlink makes the node fail with
  `File exists`, and must be removed before rerunning.

Run these with `SAGE_BENCH_MANIFEST=manifest_f9468.json` and
`docs/ai/simd_benchmark/run.py LABEL BINARY CONFIG`; the F9468 manifest and its
`f9468_bucket<N>.json` config variants live in `target/simd-benchmark`.

### Pre-SIMD confirmation

The result does not depend on any SIMD work. The frozen pre-SIMD baseline
(`42510e4`, AoS `Theoretical` records, no AVX2, `d8e5adf…`) gains **4.15x** on
F9468 and **4.24x** on F9477 from `bucket_size` alone. The commit that
introduced columnar pages and AVX2 is worth only 1.09-1.15x at the production
page size, which is why the original 11.9% measurement for it was so small.

### Caveats

Both datasets share a narrow search regime. Nothing here tests wide-window or
open search, where the precursor-filtered slice grows with tolerance and large
buckets are expected to invert. One run per configuration.

## Open follow-ups

Ranked by expected value given where time now goes. Profile before optimizing
any of them -- the one lesson this whole exercise keeps repeating is that the
obvious target was the wrong one.

### Profile scoring (largest remaining unknown)

`page_search_batch` is now only ~39% of search CPU, down from ~87%. **The other
~61% is scoring and has never been attributed.** Instrument `score_candidate`,
`build_predicted_dense` and the `matched_peaks_with_isotope` preliminary pass
before touching any of them.

One concrete candidate to confirm or dismiss with that profile:
`matched_peaks_with_isotope` allocates and zero-fills
`vec![PreScore::default(); potential]` once per isotope error, the RT/IIM
eviction loop walks all `potential` entries, `AddAssign` concatenates all four
vectors, and `trim_hits` heapifies the concatenation -- four O(potential) passes
per spectrum for a handful of non-zero entries. Measured candidate counts are
578 mean / 1,182 p95 per query (summed over the four isotope errors), so this is
real but bounded work. Cheap to fix, payoff unknown.

### Re-evaluate the mass-scan kernel in the new regime

At `bucket_size = 4,194,304` the mass scan is **58.6%** of `page_search_batch`,
up from 10.3%, and slices run ~96 fragments instead of ~0.8. The AVX2 kernel now
runs on the dominant phase with lanes that actually fill, so the old utilization
diagnostic (0.431% of ranges reaching eight masses) no longer describes it.
Re-run that diagnostic before concluding anything about the kernel -- including
the previously-skipped Scalar-SoA ablation, which is a more informative question
now than when it was proposed.

Of the four SIMD candidates proposed before any profiling existed, "SIMD across
windows on one page" is the only one this regime change revives, and only if the
re-run diagnostic supports it. The other three remain unsupported by evidence.

### Derive `bucket_size` instead of hardcoding it

4,194,304 leaves 48 pages and is near-degenerate. Returns are already flat from
2,097,152 (451 ms) to 4,194,304 (411 ms) on the F9477 sample, so an intermediate
value captures nearly all the gain with a saner structure.

Balancing bound-resolution cost `(F/B)·log2(B)` against scan cost
`peaks · B · precursor_tol` gives `B* ∝ sqrt(F / (peaks · precursor_tol))` --
so it should grow with the square root of database size and shrink with peak
count and precursor tolerance. The constants are hand-waved (plugging F9477 in
gives ~1.5e7 against a measured optimum nearer 2-4e6), so the *form* looks right
and the coefficient needs calibrating; treat it as a shape to fit, not a formula
to ship. `Builder::make_parameters` has everything needed except `F`, which is
known by the time of the `par_chunks_mut(bucket_size)` call.

**Blocking prerequisite: a wide-window/open search test.** Both benchmarked
datasets are narrow (precursor ~±6 ppm, fragment ~±13 ppm). The model predicts
large buckets *invert* as precursor tolerance grows, and Sage's index is
organized for open search in the first place. Do not change a production default
before measuring there.

### Remaining memory-parallelism headroom

Lane count is spent: 16 slots gave 1.94x on the bound phase, 32 gave 2.50x, so
the core's outstanding-miss slots are saturated. Further gains need *fewer*
misses, not more overlap -- a two-level summary index (sample every 512th
peptide ID per page into a ~1.5 MB, L3-resident side array, then scan the
narrowed window) would cut dependent DRAM probes from ~15 to ~2.

With `bucket_size` raised this is second-order: bound resolution is 41.4% of
`page_search_batch`, which is 39% of search, so making it free gains ~16% of
search. Its real value is robustness -- it would let `bucket_size` return to a
size with good mass locality, and should hold up in the wide-window regime where
large buckets are expected to fail. Do it if open search matters or if scoring
profiling turns up nothing larger.

### Not worth doing

Reorganizing the index precursor-major (bucket by peptide ID, sort by fragment
mass within bucket) was considered and rejected on measurement. It would make
the precursor range direct arithmetic and delete the bound-resolution phase
outright, but fragment-first does **5.7x less fundamental work** at these
tolerances -- 3,200 window searches against 18,258 theoretical-fragment lookups
(578 candidates x 31.6 fragments/peptide) -- and the per-query working set is
143 KB mean / 292 KB p95, not the L1-resident figure that would justify the
trade. Fragment-first was never the problem; the memory layout beneath it was.
