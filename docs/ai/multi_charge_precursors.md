# Multi-charge precursors: every listed charge, pooled (2026-10-15)

A pmsms precursor whose charge decharging couldn't settle carries several charges
(`charges` 23, 24, 34, 234; 29% of F9477's final-pass precursors). The pmsms reader
turns it into one `Precursor` per charge, all with the same m/z, as upstream's MGF reader
does for `CHARGE=2+ and 3+`. But the scorer only ever read `precursors.first()`, so such a
spectrum was searched at its first charge only: silent data loss, not a design choice.

## Semantics (user's decision: pooled)

`Scorer::candidate_charges(&query.precursors)` returns every charge listed for the first
precursor's ion (precursors with the same m/z, bit for bit), in listed order. The search
pools them as it already pooled the configured range for a spectrum without a charge:
one preliminary top-k across charges, one rank-1 PSM, whichever charge scores best, so
the FDR keeps one PSM per spectrum. Precursors with another m/z (chimeric mzML) still
don't add their charges; no charge, or `override_precursor_charge`, still means the
configured range. Rejected: one spectrum per charge (what the pipeline's MGF/mzML
exports write), which reports several rank-1 PSMs per spectrum and lets charge states
not compete.

Charges are listed as single digits (pmsms `charges`, MGF `CHARGE=`), so
`CandidateCharges` keeps listed charges in a fixed `[u8; 9]` (more is a panic) and the
configured range separately; no allocation per spectrum.

Sliced search makes one query per listed charge. A query whose spectrum as read would
search other charges too is scored on a copy whose only precursor carries its charge,
and `rank_across_charges` merges them: per-charge top-k instead of pooled, so sliced
never scores lower than monolithic (`sliced_search.md`, Charges).

## Verification

- **Unchanged where nothing is listed twice:** G23 `.d`, the short zlib mzML and MGF test
  inputs (one charge per spectrum), and G23 with `override_precursor_charge` sliced are
  identical to `56ae675`.
- **Short pmsms test cut** (442 of 2,000 precursors with several charges): single-charge
  spectra identical except `protein_groups` (grouping is global); of the multi-charge
  ones, 33 gain a PSM, 173 keep peptide and charge with only `scored_candidates`,
  `poisson` and `delta_next` changed (pooling), 156 change peptide or charge. One scores
  lower: a charge-3 preliminary candidate pushed the old winner out of the pooled top-k,
  as SAGE's pooled search for charge-less spectra already does. Sliced against
  monolithic: single-charge identical, multi-charge never lower.
- **F9477 final pass** (regression run 11's command and inputs): 890,156 single-charge
  precursors with identical search columns (only `protein_groups` of 1,814 moved);
  of 355,897 multi-charge ones, 73,073 gain a PSM and 143,394 change (81,032 to another
  charge). SAGE at 1% FDR: 93,249 → 93,658 PSMs, 20,761 → 20,595 peptides.
- **F9477 through the pipeline** (`jobs/f9477_best.toml`, reflinked copy of the
  `56ae675` nodes, so mkpmsms output is byte-identical; regression runs 10–12 before,
  16–18 after; peptidoforms / ions at 1%):
  - recalibration-pass SAGE: 10,079 / 10,712 → 10,250 / 10,926;
  - final SAGE: 20,761 / 22,701 → 20,705 / 22,717;
  - **mokapot: 28,417 / 31,027 → 29,120 / 31,778 (+2.5% / +2.4%)**, spectra
    117,752 → 121,372. Above mokapot's input sensitivity (a few hundred either way).
  - Cost: final `run_sage` node 147 → 161 s; its search including reading spectra
    (Stage 1 of streaming moved reading into the search) ~101 → 124 s, as 29% of
    precursors are searched at two or three charges. A first, direct run had been
    reported as "search 72 → 137 s", which compared a search without reading to one
    with it.

## mzML: "possible charge state" (MS:1000633)

mzML has a term for this case: MS:1000633 "possible charge state", repeated in
`selectedIon`, one per candidate ("known to be one of several possible values"), next
to MS:1000041 "charge state" for a determined one. The mzML reader now reads it: a
precursor's charges are its charge state and every possible charge state, in document
order without repeats (the set of all, when both are present), each a `Precursor` with
the same m/z, so the pooled search above applies (`mzml::one_per_charge`,
`possible_charge_states_join_the_charge_state`). Files without MS:1000633 read as before.

`pmsms2mzml --possible-charges` (default off; whether FragPipe reads it is unchecked)
writes a multi-charge precursor as one spectrum with one MS:1000633 per charge instead
of one spectrum per charge. The short F9477 test cut written that way and searched by
SAGE gives the same rank-1 peptide, charge and hyperscore as the pmsms search for all
1,862 spectra with a hit (408 multi-charge), and the same 905 peptides at 1%.
