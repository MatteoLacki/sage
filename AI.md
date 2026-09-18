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

## Code style for agents

Prefer readable, C/C++-style imperative control flow — sequential `?`-early-returns,
plain `if`/`match`, simple loops — over deeply chained combinators (`.and_then()`,
`.map()`, `.filter()`) when a sequence chains multiple independent fallible steps,
especially from different data sources. `?` and `.and_then()` chains compile to
identical code (verified for this fork's `Scorer::build_predicted_dense`,
`scoring.rs`) — there's no performance tradeoff, so pick whichever reads top-to-bottom
like a checklist, not whichever chains most tersely. Iterator adapters over a genuine
sequence (`.map()`/`.filter()`/`.fold()` over an actual collection) are fine; a
combinator chain used only to avoid writing an early-return is not.

Prefer naming a concept over commenting it. Before writing a comment that explains
what a condition means or why a check is safe, check whether extracting a
well-named function/method/field would make the comment unnecessary — the name
then carries the explanation at every call site, not just the one currently being
written, and can't drift out of sync with the logic the way a comment can. Example
(`spectrum.rs`/`scoring.rs`, 2026-09): `query.peak_charges.is_empty()` guarded by a
14-line comment explaining "this means deisotope=false, which means X, which means
Y is safe" at its one call site was replaced with `query.is_deisotoped()`, a method
whose own doc comment carries that explanation once, derived from the existing
field so it can't desync from it. Still write a comment for a genuine non-obvious
*why* (a hidden constraint, a workaround, a surprising invariant) that a name alone
can't carry — this is about replacing explanatory comments with structure, not
eliminating all comments. See the `self-documenting-naming` skill
(`~/.claude/skills/self-documenting-naming/SKILL.md`).

## Documentation index

This file stays a short overview; design rationale/history for each feature
lives in `docs/ai/`, one file per topic:

| File | Covers |
|------|--------|
| `docs/ai/pmsms_input.md` | `pmsms`/precursors binary input, per-precursor ppm tolerance, per-fragment (spline) ppm tolerance, MS1/MS2 intensity PIN columns |
| `docs/ai/predicted_rt_iim.md` | `--predicted-rt`/`--predicted-iim` hard-eviction filtering, dense peptide-index lookup, external RT/IIM as LDA features, `combined_score` soft ranking penalty |
| `docs/ai/predicted_fragment_intensity.md` | `--predicted-fragment-intensity-*`: optional MS2 fragment-intensity reader (job-scoped pointer parquet + shared `arrays.mmappet`), feature-only (no hard filter), `ms2_entropy_similarity` |
| `docs/ai/unimod.md` | `[UNIMOD:<id>]` modification notation support (`crates/sage/src/unimod.rs`) |
| `docs/ai/dump_peptides.md` | `dump_peptides` binary's mass-sorted output (for a separate consumer); necroflow's structural (not content-addressed) provenance hashing |
| `docs/ai/dump_fragment_index.md` | `dump_fragment_index` binary: ppm-binned precursor-mass x fragment-mass occupancy grid of the fragment index, mmappet flat-column output (diagnostic only) |
| `docs/ai/simd.md` | Fragment-index search performance: SoA/AVX2 page scan, interleaved peptide-ID bound resolution, `bucket_size` sizing, and the F9477/F9468 benchmarks |
| `docs/ai/reuse_index_bins.md` | `page_search_batch`: sharing one page's precursor-scoped range across a spectrum's fragment windows |
| `docs/ai/dumped_peptides_positional_predictions.md` | `--predicted-rt`/`--predicted-iim`/`--predicted-fragment-intensity-index`: positional (`peptide_row`) instead of `sequence`-keyed, `dumped_peptides_sha256` fingerprint safety net |

Check `summarise/`-style freshness only matters for the top-level monorepo;
within this vendored fork, treat each `docs/ai/*.md` file as current unless
`git log` on the files it names says otherwise.

## Search performance and `bucket_size`

Search time on narrow (non-open) searches is dominated by **memory latency in the
fragment index, not arithmetic**. Locating the precursor-compatible peptide-ID
range inside each mass page is a dependent chain of ~15 DRAM probes, and it was
94.1% of search cycles; the AVX2 mass-comparison kernel only ever touched the
other 5.9%.

Two things fix it, and the config one dominates:

| | F9477 | F9468 |
|---|---:|---:|
| `bucket_size` 32,768 -> 4,194,304 (config only) | 7.53x | 7.43x |
| Interleaved bound resolution (`resolve_page_peptide_bounds`) | 2.27x alone | 2.21x alone |
| Combined, full-dataset search | **9.37x** | **9.93x** |

Both verified output-exact (all non-ID columns and every fragment annotation).

**`bucket_size` is the single highest-leverage knob here, and upstream's default
of 8,192 (`Builder::make_parameters`) is the worst value measured.** Bigger pages
win because page-group count, not page length, sets the probe budget: fragment
windows are narrow and land on ~1 page whatever its size, so larger pages mean
fewer groups and more window reuse per group. It is a *build-time* property --
the index is mass-sorted then `par_chunks_mut(bucket_size)`-ed and each chunk
sorted by peptide ID -- so it cannot adapt at search time.

Caveats before changing a production config: both datasets used narrow
tolerances (precursor ~±6 ppm, fragment ~±13 ppm, `max_fragment_charge: 1`).
Large buckets are expected to *invert* for wide-window/open search, where the
precursor-filtered slice grows with tolerance. Untested there. Optimal size
should scale roughly as `sqrt(F / (peaks * precursor_tol))`, so it belongs
derived from index size rather than hardcoded.

Corollary worth knowing when benchmarking: several optimizations here are
*starved* by a badly-sized `bucket_size` and measure as noise at 32,768.
On F9468, index sharing (`062f7b3`) is slightly **negative** at 32,768
(738.6 -> 748.9 s) but worth 1.49x at 4M, and the columnar/AVX2 change is worth
1.15x at 32,768 but 2.07x at 4M; F9477 shows the same pattern (-2% / +22% and
1.09x / 1.94x). Always measure page-search changes at more than one
`bucket_size` -- a single measurement at 32,768 will call a real win noise.

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
