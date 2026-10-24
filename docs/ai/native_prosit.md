# Native Prosit predictions from compact peptides

`predict_prosit` streams a stored Sage compact peptide table through a native
Prosit 2023 timsTOF TensorRT engine and writes job-scoped fragment predictions
that Sage can mmap directly. This is an optional Linux executable; ordinary
Sage builds do not link TensorRT or CUDA.

## Data path

```text
FASTA + compact row columns + digestion/modification configuration
  → CompactView → Rust token IDs, charge one-hot, scaled CE
  → C ABI → TensorRT engine → Dense FP16 [N,29,6] (on GPU)
  → CUDA validity filtering, top-K, Sage annotation IDs, optional sqrt8
  → packed host records → Rust column writer
  → arrays.mmappet + positional index.mmappet
```

One preparation thread, one inference thread, and one writer overlap through
three reusable batch buffers. Memory used for batches is bounded independently
of peptide count. The existing compact table is mapped, and the proteome is
loaded once. Neither peptide strings nor a full prediction tensor are created
for the entire table. Batches remain necessary: the deployed engine accepts
1–1000 rows per call.

The GPU receives token tensors and returns only retained intensities, annotation
IDs, and counts through CUDA copies. The writer uses buffered file writes directly into mmappet columns,
then flushes and fsyncs them. It does not use writable mmap or a network protocol.
The resulting files are suitable for Sage's existing read-only mmap readers.

## Build and run

Build the Rust executable:

```sh
cargo build --release --bin predict_prosit
```

Compile the optional bridge in an environment with the matching TensorRT 8.6
and CUDA compiler (`nvcc`), headers, and libraries:

```sh
bash crates/sage-cli/native/build_prosit_bridge.sh /absolute/path/libprosit_trt.so
```

The build script accepts `NVCC`, `TENSORRT_INCLUDE_DIR`,
`TENSORRT_LIBRARY_DIR`, and `CUDA_HOME`. The deployed plan used during validation
requires TensorRT 8.6 and CUDA 12. Use the same engine and compatible GPU/runtime
as the serving deployment. The bridge loads dynamically; it is not a general
loader for arbitrary models. This executable requires bridge ABI 3; rebuild
the bridge alongside the executable. ABI 1/2 are rejected before model creation.
The direct-Dense path uses a newly built engine; the legacy FP16-output engine
remains supported with its original model-order semantics.

```sh
LD_LIBRARY_PATH=/absolute/runtime/directory \
target/release/predict_prosit \
  --fasta /absolute/proteome.fasta \
  --config /absolute/database.json \
  --compact-peptides /absolute/compact-table \
  --engine /absolute/model.plan \
  --library /absolute/libprosit_trt.so \
  --output /absolute/new-predictions \
  --min-charge 2 --max-charge 4 --collision-energy 30
```

`--top-k K` retains at most K most intense valid fragments for each peptide and
precursor charge; the default is **12**. `--top-k 0` retains all valid fragments
(the behavior used in the original benchmark below). Selection follows length,
precursor-charge and negative-intensity filtering. Equal intensities prefer the
earlier model fragment index; positive and negative zero tie. Retained fragments
are in descending intensity order for a direct-Dense engine; the legacy engine
retains model order. In FP16 storage mode, selected engine values are not renormalized.
The setting is recorded in both manifests. This caps stored fragments, not model
compute: TensorRT still evaluates its full 174-output tensor. No additional
fragment m/z or search-specific charge filter is implied by `--top-k`.

### Normalized sqrt8 output

`--intensity-encoding sqrt8_basepeak_v1` enables lossy storage with no saved
scale. After validity filtering and top-K selection, normalize each spectrum by
its maximum retained intensity, then encode in FP32 arithmetic:

```text
q = round_to_nearest_even(255 * sqrt(intensity / maximum))
decoded_relative_intensity = (float32(q) / 255)^2
```

Filtering, selection, annotation mapping, and encoding execute in one CUDA
kernel immediately after TensorRT on the same stream. For the direct-Dense
engine, the kernel reads the [N,29,6] output before LeakyReLU or flattening.
The engine omits those two operations and emits FP16 directly, avoiding a
redundant FP16-to-FP32 output conversion observed during engine inspection. FP16 mode uses the same
kernel without quantization. Rust only separates packed records into columns
and builds positional offsets; no CPU top-K or encoding fallback remains.
Each transferred record contains one count byte, K intensity slots, and K
annotation slots, with unused slots discarded by the writer. K=12 sqrt8
transfers 25 bytes per prediction instead of the previous 348-byte raw FP16
tensor. K=0 and K>174 use 174 slots. Invalid nonfinite output produces an error
count (255), which prevents library publication. The stored format is unchanged.

The ABI-3 bridge and direct-Dense engine are at
`dense-topk-study-20261022/runtime/libprosit_trt.so` and
`dense-topk-study-20261022/runtime/model.plan`. Use them for `--library` and
`--engine`. The matching CUDA/TensorRT runtime libraries remain under
`native-prosit-study-20261020/runtime/`. Old engines and bridges are preserved
for comparison. Local validation and throughput artifacts are under
`dense-topk-study-20261022/` (`accuracy.json`, `class-i-accuracy.json`,
`class-ii-accuracy.json`, and `summary.json`). Tryptic (48,924 spectra), class I
(11,916), and class II (11,610) samples matched the independent reference for
the new Dense output exactly. Compared with the old engine, selected fragment
sets changed in 1.55%, 1.29%, and 2.31% of spectra respectively; minimum decoded
cosines were 0.9495, 0.9794, and 0.9087. These checks do not establish unchanged
identification performance. Four alternating benchmark pairs of 750,018
predictions gave median pipeline times of 11.35 s before and 11.19 s after:
the apparent 1.4% throughput gain is within run-to-run variation. Full libraries
have not been regenerated with the new engine.

GPU fixture tests (NumPy and a CUDA device required):

```sh
nvcc -std=c++17 -O3 --fmad=false crates/sage-cli/native/test_prosit_pack.cu -o /tmp/test_prosit_pack
python crates/sage-cli/native/test_prosit_pack.py /tmp/test_prosit_pack
```

Zero/all-zero inputs produce zero codes. An empty spectrum remains empty.
Negative values are filtered and nonfinite values rejected before encoding.
Retained zero codes remain present channels, distinct from missing predictions.
The maximum is discarded: original amplitudes cannot be recovered. This is
intentionally different from the earlier proposed `sqrt8_scaled_v1` codec,
which stored an FP32 scale per spectrum. The user explicitly chose normalized,
lossy results on 2026-10-21.

The output's `arrays.mmappet/schema.txt` becomes:

```text
uint8 predicted_intensity
uint8 annotation_id
```

Offsets and annotation IDs are unchanged. `arrays.mmappet/meta.json` explicitly
declares version 2, `intensity_encoding: sqrt8_basepeak_v1`,
`normalization: per_spectrum_base_peak`, and `rounding: nearest_even`.
The same metadata is in the root/index manifests. The updated Sage reader
requires this contract for UINT8 values and decodes only the requested
fragments during scoring. Both regular and sliced searches use it without
expanding the library into RAM. Legacy FP16 files remain readable, and `fp16`
remains the executable's compatibility default; select sqrt8 explicitly.

For a full 12-fragment prediction, columns shrink from 52 bytes (24 intensity,
12 annotation, 16 offsets) to 40 bytes (12 intensity, 12 annotation, 16 offsets).
The completed tryptic library, including missing slots, shrank from
1,488,919,189 to 1,160,643,566 column bytes: 22.05%. A separate conversion copy
preserves the source data. Across 27,415,115 predicted spectra, cosine against
the normalized FP16 source averaged 0.9999970095 (worst 0.9999906485). No claim
about identification performance follows from this codec measurement alone.
The GPU still computes all 174 outputs, so this is primarily a storage change.

Workspace artifacts and the independent NumPy converter are in
`sqrt8-study-20261021/`. The large library queue remains stopped; old outputs
are not silently reinterpreted, and no live Triton server was changed.

`--config` takes the database section, as `dump_peptides --config` does. Use the
same FASTA and configuration that generated the stored compact table. A full
export stores its checked compact fingerprint as `dumped_peptides_sha256`.

Use the output in Sage with:

```text
--compact-peptides /absolute/compact-table
--predicted-fragment-intensity-index /absolute/new-predictions/index.mmappet
--predicted-fragment-intensity-cache /absolute/new-predictions/arrays.mmappet
```

The prediction charge range must match the search's configured precursor charge
range. Existing Sage positional indexing assumes this; fingerprint equality
alone does not validate charges, collision energy, or model choice.

## Model and modification semantics

Inputs are INT32 `[N,30]` tokens, FP32 `[N,6]` one-hot charges, and FP32 `[N,1]`
collision energy multiplied by `0.01`. The direct-Dense engine outputs FP16
`[N,29,6]`; the legacy engine outputs post-LeakyReLU FP16 `[N,174]`.

The tokenizer follows the deployed compact Prosit 2023 vocabulary:

- Unmodified canonical amino acids use IDs 1–20; padding is 0.
- `M[UNIMOD:35]` uses ID 21.
- `C[UNIMOD:4]`, `K[UNIMOD:259]`, and `R[UNIMOD:267]` use their base residue IDs.
- Other residue modifications, terminal modifications, noncanonical residues,
  empty peptides, and lengths above 30 are unsupported in `original` mode.

Configuration-resolved Unimod IDs establish exact mass-bit mappings. No fuzzy
mass matching is performed. Prefer explicit `UNIMOD:<id>` configuration entries.
`CompactView` supplies Sage's modification placement, protein-terminal rules,
and decoy reversal semantics.

Unsupported peptides fail the run by default. `--unsupported skip` instead
writes `-1/-1` positional entries for all charges of those peptides and records
counts in the manifest. Explicit `--modifications all` strips modifications
before encoding, matching a different prediction policy; it must be chosen
deliberately.

Prosit orders fragments by ordinal, interleaving y and b, with charges 1–3.
Sage IDs use b and y blocks indexed by N-terminal backbone position:

```text
b ordinal o, charge z:             (o - 1) * 3 + z - 1
y ordinal o, charge z, length L: 87 + (L - 1 - o) * 3 + z - 1
```

Keep ordinals below peptide length, fragment charges at most precursor charge,
and nonnegative intensities. Reject nonfinite outputs. FP16 mode preserves
signed zero and raw intensity bits; sqrt8 mode normalizes and quantizes as above.

## Output and failure behavior

Output contains:

- `arrays.mmappet`: FP16 or normalized UINT8 `predicted_intensity`, uint8
  `annotation_id`, with encoding metadata.
- `index.mmappet`: int64 `start`, int64 `end`, with dense slot order
  `peptide_row * charge_span + charge - min_charge`.
- `manifest.json` and index metadata: engine/bridge/config provenance, compact
  fingerprint, row and charge ranges, prediction policy, and stage timings.

A new output directory is required. Files are built in a sibling
`.<name>.partial-<pid>` directory and published by directory rename only after
all files and metadata have been flushed. Failures leave that partial directory
for inspection; they do not publish it. An entirely empty fragment dataset is
not published because the existing intensity reader cannot mmap empty arrays.
Failed runs currently restart from the beginning.

`--row-start` and `--row-count` are validation conveniences. Partial exports
deliberately omit Sage's full-table fingerprint and cannot be used as ordinary
full-table prediction indexes. `--reference-tsv` optionally materializes strings
and tokens for validation outside the pipeline timer.

Existing predictions are not reused or deduplicated. Reuse is deferred. The
manifest records the model and policies needed to design a future key; peptide
sequence alone would not be sufficient.

## Validation, 2026-10-20

A deterministic 20-protein human FASTA subset yielded 250,728 compact rows
(125,313 decoys; 94,788 modified rows), lengths 8–15. All Rust tokens and
supported/unsupported decisions matched the deployed Python parser.
722 terminal-modified rows were skipped in original mode. Charges 2–4 yielded
750,018 actual predictions and 30,076,081 retained fragments.

Three alternating-order runs compared native inference with the live packed
Triton endpoint at concurrency 16, using identical input order, batches, model
plan, CE, output columns, and fsync. Every output column was byte-identical in
every run.

| Median, includes output writes and fsync | Native | Packed gRPC |
|---|---:|---:|
| Pipeline seconds | 9.558 | 11.163 |
| Actual predictions/second | 78,471 | 67,185 |

This is 16.8% higher throughput (14.4% less pipeline wall time). The gRPC
baseline starts from pre-materialized strings; their export is excluded.
Native compact input preparation is included. Native process startup and model
loading are excluded from the pipeline figures; median total process time was
10.330 seconds versus the already-running gRPC service.

Median native stage work: preparation 0.059 s, inference including CUDA copies
9.488 s, packing/writing 0.696 s, final flush 0.065 s. Stages overlap, so their
times must not be added. These results do not support a 75% wall-time reduction.
They are a short-peptide fixture result on the local RTX A4000, not a 600M-row
capacity or sustained-I/O test and not a benchmark of cache reuse.

The workspace study directory `native-prosit-study-20261020/` contains the
fixture generator, comparison script, full per-round JSON, hashes, and report.
Core tests cover token vocabulary and materialized/compact equivalence;
pipeline tests cover mmap reader compatibility, fragment mapping, signed zero,
unsupported positional slots, nonfinite output, and error propagation.

## Building a direct-Dense engine

`native/export_dense_prosit.py` (under `crates/sage-cli/`) reads the existing
FP16-output ONNX export, removes the trailing LeakyReLU, flatten reshape, and
output cast, and exposes `dense_intensities` with logical shape [N,29,6]. It
leaves all trained initializers unchanged and records source/output hashes.
It requires the ONNX Python package.

`native/build_dense_prosit.py` builds with TensorRT 8.6 on the target GPU,
using FP16 tactics and the original engine's input profiles. Its default
`--output-dtype fp16` exposes the engine's natural Dense precision directly.
The optional FP32 output was inspected and found to add a redundant output
conversion in the tested engine. Do not infer every internal layer's precision
from this output format; TensorRT fuses the network.

```sh
python crates/sage-cli/native/export_dense_prosit.py original_fp16out.onnx dense.onnx
python crates/sage-cli/native/build_dense_prosit.py dense.onnx model.plan --reference-engine original.plan
```

Dense outputs are selected and sorted by raw intensity, with ties broken by
model index, before sqrt8 encoding. Metadata records `model_output`,
`fragment_order: descending_intensity`, and `native_bridge_abi: 3`.
The disk codec and Sage annotation encoding are unchanged. The legacy engine
retains `fragment_order: model`. A rebuilt engine can change low-order numeric
results and which nearly tied fragments are selected; compare by annotation ID,
not by column position, when validating against older libraries.
