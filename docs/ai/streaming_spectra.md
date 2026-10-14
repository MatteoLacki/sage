# Streaming spectra: read, process and score one spectrum at a time (2026-10-14)

Plan: `plans/streaming_spectra.md` (necromerge2 root). SAGE used to read a batch of
input files whole, keep every processed spectrum in one `Vec<ProcessedSpectrum>`, and
only then search it. Every result of scoring lives in the `Feature`, so the search
needs one spectrum at a time; holding all of them is the reader's habit, not a need.

## Stage 0: what the spectra hold

`runner.rs` logs `- memory, <checkpoint>: anonymous …, peak RSS …` (`RssAnon` and
`VmHWM`, `process_memory::log_memory`) when the database is ready, when spectra are
read, after the search and after the spectra are dropped. The batch reader also logs
the exact bytes its processed spectra own (`ProcessedSpectrum::owned_bytes`). Measured
with `56ae675` plus these logs:

| input | spectra | processed spectra own | anonymous memory while reading | peak anonymous |
|---|---|---|---|---|
| F9477 final pass (pmsms, regression run 11) | 1,246,053 | 8.62 GB | 4.74 → 12.73 GB | 14.92 GB |
| G23 `.d` | 35,351 | 0.04 GB | ~0 | 5.59 GB |
| long MGF (`temp/streaming_test_data`) | 137,411 | 1.00 GB | 4.18 → 8.11 GB | 12.12 GB |
| long zlib mzML | 137,411 | 0.98 GB | 4.29 → 8.78 GB | 8.75 GB |

- On F9477 the processed spectra are 8.6 of the 14.9 GB anonymous peak.
- MGF and mzML cost 4–8 times the spectra they produce while reading: MGF reads the
  whole file into one `String` (4.45 GB), mzML keeps every raw `f32` array of the batch
  until processing.
- Dropping the spectra after the search frees them, but anonymous memory barely falls
  (14.96 GB after dropping 8.6 GB on F9477): glibc keeps freed heap mapped. Only not
  allocating lowers the peak.

## Stage 1: pmsms and TDF read by index

`PmsmsSpectra` (`sage-cloudpath/src/pmsms.rs`) and `TdfSpectra` (`tdf.rs`) open an
input once and read spectrum `i` on demand: pmsms from the precursor table and the
mapped fragment columns, TDF through `SpectrumReader::get`. The batch readers
(`pmsms::parse*`, `TdfReader::parse`) became parallel maps over the same `read`, so
both paths build each spectrum with one code.

`Runner::process_chunk` searches a batch with `search_chunk_by_index` when every file in
it can be read by index (`IndexedSpectra`):

    (0..n).into_par_iter().map(|i| score(process(read(i)))).collect()

rayon walks each thread's index range lazily, so about one spectrum per thread is
alive, and `collect` keeps index order: features come out in today's order. Files of
a batch are searched one after another and appended in file order, as before.

The batch reader stays for:
- **TMT and LFQ runs:** TMT quantifies from the spectra after the search, LFQ needs
  the MS1 spectra (Stages 2 and 4).
- **mzML and MGF:** Stage 2 gives them offset indices.
- **Sliced search and prefilter:** sliced search moves in Stage 3; prefilter stays.

A TDF file that fails to open is logged and contributes nothing, as in the batch
reader; a pmsms read error panics, as before.

### Verification (Stage 1 against `56ae675`)

Compared PSM by PSM (every column except `psm_id` and `sage_discriminant_score`) plus
all matched fragments; scripts in `temp/streaming/` (necromerge2 root):

| input | path | result | peak anonymous memory |
|---|---|---|---|
| F9477 recalibration pass (regression run 10) | by index | identical (217,366 PSMs) | 4.98 GB |
| F9477 final pass (regression run 11) | by index | identical (733,681 PSMs, 4,284,005 fragments) | 6.45 GB, was 14.92 GB; wall 137 s, unchanged |
| G23 `.d` | by index | identical (16,661 PSMs) | 5.62 GB, unchanged (its spectra are 0.04 GB) |
| `temp/streaming_test_data` short and long: pmsms | by index | identical | |
| short and long: zlib and uncompressed mzML, plain and `.gz`; short MGF | batch | identical | |
| long MGF, plain and `.gz` | batch | `posterior_error` and q-values differ | |
| numpress mzML | batch | crashes, as before (Stage 2) | |

The long MGF difference is SAGE's own run-to-run noise: `56ae675` run on the plain and
the gzipped file (identical spectra) differs in the same columns
(`sage_discriminant_score`, `posterior_error`, every q-value); scores, features and
fragments are identical. On F9477's pmsms the same LDA noise never reaches the q-values.
