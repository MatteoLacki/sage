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

## Stage 2a: MS-Numpress decoding

The pipeline's only mzML flavour (`pmsms2mzml --numpress`) crashed SAGE: the numpress
cvParam was unknown, the bytes were read as float64, and `copy_from_slice` panicked on
the remainder. `sage-cloudpath/src/numpress.rs` decodes the three codecs in safe Rust:
a bounds-checked half-byte reader, `decode_linear` (m/z), `decode_pic` and
`decode_slof` (intensities), appending f32 to a reused buffer; corrupt or truncated
input is an error, never a panic (the linear predictor's i64 arithmetic is checked
too: the reference overflows silently).

No external crate: `numpress` 1.1.0 has the linear codec only; `numpress-rs` 1.1.1 (what
`mzdata` uses, via `mzdata-bindata`) is raw-pointer unsafe code that panics on arrays
under 8 bytes, over-reserves `decode_pic` 40 bytes per input byte and returns `Ok(0)`
from it, and `shrink_to_fit`s every call.

The mzML reader maps MS:1002312/1002746 to linear, MS:1002313/1002747 to pic,
MS:1002314/1002748 to slof; the 10027xx terms ("followed by zlib") also set the zlib
flag, so zlib is inflated once even where it is also declared as MS:1000574 (as
`pmsms2mzml --numpress` does). The codec resets at each `binaryDataArray` and
overrides the float-width term.

Verified:
- **Reference vectors:** arrays encoded and decoded by the reference C++ (vendored in
  `git/pmsms2mzml`) decode to the same f32 bits; every truncated prefix errors or decodes
  a prefix of the values; random bytes never panic.
- **Real data** (`numpress_mzml_matches_its_zlib_twin`, ignored test, `SAGE_NUMPRESS_PAIR`):
  the numpress and zlib mzML of `temp/streaming_test_data` decode to identical arrays,
  intensities and m/z alike (short: 2,498 spectra, 6.6M peaks; long: 137,411 spectra,
  367.6M peaks; m/z difference 0 ppm: pmsms2mzml's optimal fixed point round-trips every
  f32).

## Stage 2b: one-spectrum parsers

- **mzML:** the parser's state lives in a `SpectrumBuilder` fed XML events. The async
  whole-document `MzMLReader::parse` (still used for LFQ, TMT and remote files) and the
  new synchronous `parse_spectrum(&[u8])` (one `<spectrum>` element; text after its
  `</spectrum>` ignored) drive the same builder. zlib is inflated with synchronous
  `flate2` (already in the build through `async-compression`).
- **MGF:** `MgfReader::parse_header` (the lines before the first `BEGIN IONS`) and
  `parse_block` (one `BEGIN IONS` … `END IONS`); `parse` is header plus blocks.

Both used to leak state from one spectrum into the next; per-spectrum parsing can't
reproduce that, so every spectrum now starts clean in both paths:
- **mzML:** isolation-window offsets, the compression flag, the float width and an
  unpushed precursor (selected-ion m/z 0) carried over. A spectrum without an isolation
  window inherited the previous one's.
- **MGF:** the file-level `CHARGE=`/`TOL=`/`TOLU=` never reached the first spectrum (its
  query state was not initialised from them), and lines between `END IONS` and the next
  `BEGIN IONS` went into the next spectrum.

Unit tests pin both (`spectra_parse_alone_as_in_their_file`,
`blocks_parse_alone_as_in_their_file`). SAGE ignores mzML `referenceableParamGroup`s, so
mzML needs no file-level context.

## Stage 2c: mzML and MGF read by index

`text_spectra::TextSpectra` opens a local mzML or MGF file and records where each
spectrum's text starts (`<spectrum` followed by whitespace or `>`; `BEGIN IONS` first on
its line). Spectrum `i` is the text from start `i` to start `i + 1` (or the end),
parsed by the one-spectrum parsers. mzML is read with an MS-level filter of 2, so MS1
spectra are skipped before their arrays are decoded.

- **Plain files:** memory-mapped; the starts are found by a parallel `memchr::memmem`
  search over 8 MB ranges.
- **Gzipped files** (`gzip_index.rs`): one sequential decompression pass finds the starts
  in the stream and records access points every 4 MB of output, checking every member's
  CRC-32 and length. zlib's zran method resumes with `inflatePrime` and
  `inflateSetDictionary`; the build has no C zlib, only `flate2`'s backend
  `miniz_oxide`, whose `DecompressorOxide` is `Clone`: an access point is a clone of it
  (bit buffer included) plus its 32 KB output ring, ~43 KB. A read resumes from the
  last point before its offset; each rayon task (`map_init`) keeps a cursor that just
  inflates forward to its next spectrum. Multi-member files (`pigz`, `bgzip`) are
  followed across members.
- **Remote URLs** (S3, Azure) keep the batch reader, which streams them.

`memchr`, `crc32fast` and `miniz_oxide` were already in the build; they are now direct
dependencies of `sage-cloudpath`, with no new crate.

The cursor keeps its inflater between reads and serves a read's start from the 32 KB
output ring when it lies there: a read inflates in chunks of up to 32 KB, so it runs past
its end, and the next spectrum starts inside those bytes. Without that, every spectrum
restarted from an access point (long gzipped MGF: 100 s of search instead of 19 s);
`consecutive_reads_continue_without_jumping` pins it.

### Verification (Stage 2 against `56ae675`)

Same comparison as Stage 1; scripts `temp/streaming/verify_stage.sh`,
`verify_stage2.sh`, `verify_gzip.sh` (necromerge2 root):
- **F9477 both passes, G23 `.d`, every pmsms input:** identical.
- **mzML zlib and uncompressed, plain and `.gz`; short MGF, plain and `.gz`:** identical,
  now read by index.
- **Long MGF, plain and `.gz`:** only `posterior_error` and the q-values differ, SAGE's
  own run-to-run LDA noise on this input (see Stage 1).
- **Numpress mzML, plain and `.gz`, short and long:** newly readable; identical to
  `56ae675`'s search of the zlib mzML of the same spectra.
- **Whole-file parsers** (LFQ run, short zlib mzML and MGF): identical.

Long test inputs (137,411 spectra), `56ae675` against Stage 2:

| input | peak anonymous | wall (search only, when measured) |
|---|---|---|
| MGF | 12.12 → 5.33 GB | 106.8 → 29.1 s |
| zlib mzML | 8.75 → 5.17 GB | 69.0 → 26.7 s |
| MGF `.gz` | | 80 → 47 s (search 18.6 s; plain 14.3 s) |
| zlib mzML `.gz` | | 73 → 36 s |
| uncompressed mzML `.gz` | | 73 → 43 s |
| numpress mzML `.gz` | | crashed → 36 s |

The `.gz` walls are from runs that also annotate fragments and write the PIN (the
`sage_checks.sh` settings), on both sides. What's left of the anonymous peak is the
database and predictions (~4.2 GB) plus the results.
