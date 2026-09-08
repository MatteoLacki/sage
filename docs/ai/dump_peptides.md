# `dump_peptides` sorts output by monoisotopic mass (2026-09-07)

`crates/sage-cli/src/bin/dump_peptides.rs` digests a FASTA into every
target/decoy `Peptide` and writes them to parquet
(`sage_cloudpath::parquet::serialize_peptides`) — see that binary's own
top-of-file doc comment for what it's for (feeding `git/featureprediction`'s
RT/IIM/fragment-intensity prediction, no search/fragment-index involved).

Added one line right after `params.digest(&fasta)`:
`peptides.sort_by(|a, b| a.monoisotopic.total_cmp(&b.monoisotopic));` —
ascending by `Peptide.monoisotopic` (`f32::total_cmp`, NaN-safe, stable).
For a **separate, SAGE-unrelated** consumer that wants the peptide table
mass-sorted; nothing internal to this binary or to `sage` itself depends on
digest order, so this is a self-contained, zero-risk change.

## Downstream: necroflow provenance is structural, not content-addressed

Learned the hard way while doing this: rebuilding this binary does **not**
give `dump_peptides`'s necroflow node a new provenance hash. Node identity
is `f(rule identity, parent *node keys*)`, not parent content — a parent's
node key (here, `source_dump_peptides_binary/<hash>/dump_peptides`, a
symlink into `target/release/dump_peptides`) stays the same string across a
plain rebuild at the same path, so the *content*-hash check
(`consumed_sha256` vs current) instead marks the child `STALE` for an
**in-place** rerun at its existing path — not a new sibling directory (see
`necroflow explain --json` on the actual job; `parent_content_changed`
reason). Confirmed empirically for `predict_fragment_intensity` too: after
regenerating `peptides.parquet` in place with the newly-sorted binary,
rerunning that node's real fill command against the *same* cache directory
made **zero new Koina calls** (row count unchanged, 714,530,556) — the
peptide *set* is unaffected by this change, only row order, and the cache
is keyed by sequence, not row position. See `git/featureprediction`'s
`docs/ai/fragment_cache_mz_sort.md` for the full verification.
