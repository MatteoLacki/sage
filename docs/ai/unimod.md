# Unimod modification support (`crates/sage/src/unimod.rs`)

Found (2026-08-21) empirically against the live Koina server: both
`Chronologer_RT` and `IM2Deep` silently ignore this fork's own `[+mass]`
mass-delta bracket notation for modifications — `PEPTC[+57.0216]DEK` and
`PEPTCDEK` come back with byte-identical predictions, only genuine
`[UNIMOD:<id>]` notation actually changes anything. Traced to Koina's own
`IM2Deep_Preprocess_AC`/`Deeplc_Preprocess_onehot` server-side code
(`wilhelm-lab/koina` GitHub): internal mods are silently regex-skipped
(no error), N-terminal mods used to crash outright (`KeyError: '+42'` —
this fork's own N-term notation, `[+mass]-SEQUENCE` with a hyphen, *does*
get split correctly as a terminal-mod segment, then the raw mass string
gets looked up as if it were itself a literal Unimod ID).

`static_mods`/`variable_mods` (`database` config) now accept a modification
value as either a plain float (unchanged) or a `"UNIMOD:<id>"` string,
resolved to that id's exact monoisotopic mass at config-load time
(`crates/sage-cli/src/input.rs`'s `resolve_unimod_refs`/
`resolve_unimod_refs_in_database` — walks the raw config JSON *before*
deserializing into `Builder`, so `Builder.static_mods: HashMap<String, f32>`
itself never changes shape; the rest of the engine only ever sees plain
floats, exactly as before). A plain float that coincides (within 5 mDa —
`unimod::COINCIDENCE_TOLERANCE_DA`) with a real Unimod entry is a config
**error**, not silently accepted — e.g. `57.0216` (what real job configs
actually write for Carbamidomethyl) errors naming `UNIMOD:4`, since letting
it through would leave that modification permanently unable to round-trip
as `[UNIMOD:4]` on output despite plainly being that modification.

Resolution is against an *active* table — `crates/sage/data/unimod.csv`
(`id,name,mono_mass`, regenerate via `scripts/build_unimod_table.py`,
same source `https://www.unimod.org/obo/unimod.obo` Koina's own
preprocessing downloads) embedded into the binary via `include_str!` and
used by default, or `--unimod-db-path <path>` (both the `sage` binary and
`dump_peptides`) to override it entirely with an external CSV of the same
shape. Either way, parsed once, lazily, into a process-global
`OnceLock` — correct for the real one-shot-CLI-process use case; test
code must take care not to have more than one test per binary actually
resolve a `UNIMOD:<id>` reference or override the active table (same
constraint, and same reasoning, as the existing capturing-logger tests in
`crates/sage-cli/src/input.rs`).

`Peptide::Display` (`crates/sage/src/peptide.rs`) round-trips: a
modification whose mass is bit-identical to one resolved via `UNIMOD:<id>`
*this run* prints `[UNIMOD:<id>]` instead of `[+mass]` — deliberately
provenance-scoped to just the ids this run's config actually referenced
(`unimod::set_reverse_table`/`lookup_reverse`), not a full ~1560-entry
reverse lookup, so this never fires on a coincidental match the config
didn't explicitly ask for. `dump_peptides`' own output, `results.sage.tsv`'s
`peptide` column, and this fork's own `--predicted-properties` `HashMap`
lookup key (`scoring.rs`) all stay consistent with each other automatically
as a result — no call-site changes needed anywhere `Peptide::to_string()`
was already used. When no `UNIMOD:<id>` reference is used at all, output is
byte-for-byte identical to before this existed.

Verified against the real, live Koina server that exposed the bug (not
just unit tests): real `dump_peptides` output for
`[UNIMOD:1]-MPEPTC[UNIMOD:4]DEK` (N-terminal + internal mod, the exact
shape that used to crash) round-trips through both `Chronologer_RT` and
`IM2Deep` without error, and the modified/unmodified predictions now
genuinely differ (e.g. IM2Deep CCS 338.67 → 346.43 Å² for the internal
mod alone).
