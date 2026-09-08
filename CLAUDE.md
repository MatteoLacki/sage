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

## Documentation index

This file stays a short overview; design rationale/history for each feature
lives in `docs/ai/`, one file per topic:

| File | Covers |
|------|--------|
| `docs/ai/pmsms_input.md` | `pmsms`/precursors binary input, per-precursor ppm tolerance, per-fragment (spline) ppm tolerance |
| `docs/ai/predicted_rt_iim.md` | `--predicted-rt`/`--predicted-iim` hard-eviction filtering, dense peptide-index lookup, external RT/IIM as LDA features, `combined_score` soft ranking penalty |
| `docs/ai/predicted_fragment_intensity.md` | `--predicted-fragment-intensity-*`: optional MS2 fragment-intensity reader (job-scoped pointer parquet + shared `arrays.mmappet`), feature-only (no hard filter), `ms2_entropy_similarity` |
| `docs/ai/unimod.md` | `[UNIMOD:<id>]` modification notation support (`crates/sage/src/unimod.rs`) |
| `docs/ai/dump_peptides.md` | `dump_peptides` binary's mass-sorted output (for a separate consumer); necroflow's structural (not content-addressed) provenance hashing |

Check `summarise/`-style freshness only matters for the top-level monorepo;
within this vendored fork, treat each `docs/ai/*.md` file as current unless
`git log` on the files it names says otherwise.

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
