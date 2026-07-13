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
- Native support for reading ionmaiden's `pmsms`/`tof2mz`/precursors binary
  formats directly, instead of only mzML/MGF/(Bruker) TDF.

## pmsms/tof2mz/precursors input

Core reader: `crates/sage-cloudpath/src/pmsms.rs`. Three inputs, given as
explicit paths (never assumed to live together in one fixed-filename
directory):

- `pmsms.mmappet/` — columnar binary, `0.bin`=tof(u32), `1.bin`=intensity(u32)
- `tof2mz.mmappet/` — `0.bin`: f32 array mapping tof_index → m/z
- precursors, either `.parquet` (via the `parquet` crate) or `.mmappet/`
  (schema.txt + one `<idx>.bin` file per column, columns looked up by name —
  see `git/mmappet`'s Python writer for the on-disk format). Needs 7 columns:
  `precursor_idx`(u64), `mz`(f64), `rt`(f64), `inv_ion_mobility`(f64),
  `charges`(i64), `fragment_spectrum_start`(u64), `fragment_event_cnt`(u64).

Two ways to pass these three paths to the `sage` CLI:

1. **`--pmsms <path> --tof2mz <path> --precursors <path>`** (preferred). Must
   be given all three together — validated in
   `crates/sage-cli/src/input.rs::Input::from_arguments`. Used instead of the
   positional spectra-paths argument; `PmsmsPaths`
   (`crates/sage-cloudpath/src/util.rs`) carries the three paths through
   `Search`, and the runner's per-file loop
   (`crates/sage-cli/src/runner.rs::process_chunk`) branches on
   `parameters.pmsms_paths` to call `read_pmsms_explicit` directly, bypassing
   the normal `FileFormat`-suffix dispatch.
2. **Positional `<dir>.pmsms`** (legacy, still works, unchanged behavior) — a
   directory whose name ends in `.pmsms`, containing `pmsms.mmappet/`,
   `tof2mz.mmappet/`, and `precursors.parquet` at fixed names. Detected by
   `crates/sage-cloudpath/src/util.rs::FileFormat::from` and read via
   `read_pmsms`, which just joins the three fixed names and calls the same
   `pmsms::parse`.

`necromerge2`'s `run_sage` Snakemake-via-necroflow rule
(`git/ionmaidentools/src/ionmaidentools/pipelines.py`) uses option 1 — no
staging directory, no symlinks, no mmappet→parquet conversion needed anymore.

## Test fixture

`crates/sage-cloudpath/tests/data/pmsms_fixture/` is a small (~400K, 10-row)
real pmsms/tof2mz/precursors triplet (both `.parquet` and `.mmappet`
precursors versions) used by `pmsms.rs`'s unit tests to assert both precursor
formats parse to identical spectra. **This repo's `.gitignore` has blanket
`data/` and `*.txt` rules** that would otherwise silently drop this fixture
(the `.gitignore` has an explicit negation carve-out for
`crates/sage-cloudpath/tests/data/` — if you add fixtures elsewhere under a
directory named `data` or with `.txt`/`.json`/`.csv`/`.tsv` files, check
`git status`/`git check-ignore -v` before assuming they're tracked).

## Testing

```
cargo test -p sage-cloudpath --lib                     # default features (no parquet)
cargo test -p sage-cloudpath --lib --features parquet  # matches the real `sage` binary build
```

`sage-cli` always builds `sage-cloudpath` with `features = ["parquet"]`
(see `crates/sage-cli/Cargo.toml`), so the actual `sage` binary always has
pmsms support; the `parquet`-feature-off path only matters if `sage-cloudpath`
is ever used as a library without it.
