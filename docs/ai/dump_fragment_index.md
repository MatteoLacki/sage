# `dump_fragment_index`: ppm-binned occupancy of the fragment index (2026-09-16)

`crates/sage-cli/src/bin/dump_fragment_index.rs` builds a real
`IndexedDatabase` and dumps a 2D histogram of its (peptide monoisotopic
mass, fragment neutral mass) pairs, for visualising where the index
actually puts mass. It is a diagnostic/plotting tool: nothing in the search
path reads its output. `git/ionmaidentools` runs it as the request-only
`bin_fragment_index` rule and renders it with necromerge2's
`scripts/mass_grid_heatmap.py` (the same renderer serves the
Koina-prediction grid, so `grid.json` carries the axis labels) -- see that
repo's `docs/ai/mass_grid_plots.md`.

## Why the histogram is accumulated in Rust, not dumped as pairs

The human FASTA at the pipeline's own `dump_peptides` settings gives 6.3 M
peptides / 199 M fragments (99.5 M pairs once decoys are dropped). Dumping
raw pairs would be ~1 GB of f32s that the consumer then has to re-bin, so
the binary bins and counts in one pass instead and writes only the grid.

## Why the grid is folded, not a true 5 ppm x 5 ppm array

Both axes are geometric: bin `k` covers `[m0 * r^k, m0 * r^(k+1))` for
`r = 1 + ppm * 1e-6` and `m0` the integer floor of the axis minimum, so a
bin is a constant ppm wide and bin index is linear in `log(mass)` (which is
what makes the rendered image a log-log plot for free). At 5 ppm that is
379 425 precursor x 627 607 fragment bins = 2.4e11 cells — 238 GB as
`u32`, and `u8` would overflow. So `--target-pixels` (default 8192) picks
the smallest **integer** fold per axis (47 and 77 for the numbers above),
which keeps every pixel edge on an exact ppm-bin boundary; `grid.json`
records `m0`/`fold`/`n_bins` so the mass axes are reconstructible.

An exact-5 ppm *sparse* triplet table was the alternative, at ~0.9 GB for
~99 M occupied cells; the folded grid is 263 MB and no zoom level a raster
can show resolves 5 ppm anyway.

## Output format

An mmappet directory in the single-flat-column convention `timstofu`'s
`binary/array_serialization.py` reads: `schema.txt` (`uint32 x`), `0.bin`
(C-order counts, precursor-major), `shape.txt` (`<nx>x<ny>`), plus a
`grid.json` sidecar for the bin parameters. mmappet has no N-D or ragged
column type and no metadata slot, hence the two sidecars — same shape as
the pipeline's own `sample_tensors.mmappet`.

## What it reads from the config

The same `database` section `dump_peptides` takes, plus the keys the
fragment set itself depends on (`ion_kinds`, `min_ion_index`).
`bucket_size` is accepted but irrelevant: it only reorders fragments within
an index page, and the histogram is order-independent, so the pipeline
leaves it out of this rule's config. Note `fragment_min_mz`/`fragment_max_mz` do **not** apply:
they filter observed peaks in `SpectrumProcessor`, so the index holds
fragments outside them (171–3943 Da for the settings above, against a
configured 200–1700 window). Fragment masses are neutral and
charge-agnostic — one entry per (peptide, kind, ion index), per
`build_from_peptides`' "all theoretical fragments are
monoisotopic/uncharged".
