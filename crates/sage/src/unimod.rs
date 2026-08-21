//! Unimod modification reference table: `UNIMOD:<id>` -> monoisotopic mass
//! delta, and the reverse for round-tripping on output. See
//! `necromerge2`'s `plans/`: Koina's RT/IIM predictors (Chronologer,
//! IM2Deep) silently ignore SAGE's own `[+mass]` bracket notation for
//! modifications -- only genuine `[UNIMOD:<id>]` notation is recognized.
//! This module lets `static_mods`/`variable_mods` reference a modification
//! as `"UNIMOD:<id>"` (resolved to the table's exact float at config-load
//! time -- the rest of the engine only ever sees plain `f32`, unchanged)
//! and lets `Peptide`'s `Display` impl print `[UNIMOD:<id>]` back for any
//! modification that *was* resolved that way this run.

use std::collections::HashMap;
use std::sync::OnceLock;

/// Same source Koina's own `IM2Deep_Preprocess_AC/1/modifications.py`'s
/// `Unimod` class downloads and parses (`https://www.unimod.org/obo/
/// unimod.obo`) -- using the identical source guarantees our resolved
/// masses agree with what Koina itself resolves for the same
/// `UNIMOD:<id>` reference. Regenerate via `scripts/build_unimod_table.py`.
const EMBEDDED_UNIMOD_CSV: &str = include_str!("../data/unimod.csv");

/// `id,name,mono_mass` -> `id -> mono_mass`. `name` is parsed and
/// discarded; kept in the CSV purely for human readability/diffing.
fn parse_csv(csv: &str) -> Result<HashMap<u32, f32>, String> {
    let mut table = HashMap::new();
    for (line_no, line) in csv.lines().enumerate().skip(1) {
        if line.is_empty() {
            continue;
        }
        let mut parts = line.splitn(3, ',');
        let id = parts
            .next()
            .ok_or_else(|| format!("unimod csv line {line_no}: missing id field"))?;
        let _name = parts.next();
        let mass = parts
            .next()
            .ok_or_else(|| format!("unimod csv line {line_no}: missing mono_mass field"))?;
        let id: u32 = id
            .parse()
            .map_err(|_| format!("unimod csv line {line_no}: invalid id `{id}`"))?;
        let mass: f32 = mass
            .trim()
            .parse()
            .map_err(|_| format!("unimod csv line {line_no}: invalid mono_mass `{mass}`"))?;
        table.insert(id, mass);
    }
    Ok(table)
}

fn embedded_table() -> &'static HashMap<u32, f32> {
    static TABLE: OnceLock<HashMap<u32, f32>> = OnceLock::new();
    TABLE.get_or_init(|| {
        parse_csv(EMBEDDED_UNIMOD_CSV).expect("embedded crates/sage/data/unimod.csv is malformed")
    })
}

/// The forward table actually in use this run: the embedded default,
/// unless `set_active_table` was called with an override (`--unimod-db-path`).
static ACTIVE_TABLE: OnceLock<HashMap<u32, f32>> = OnceLock::new();

/// Load `path` (a Unimod CSV, same `id,name,mono_mass` shape as the
/// embedded default) and make it the active table for this run, in place
/// of the embedded default. Must be called (if at all) before any
/// `UNIMOD:<id>` reference is resolved. Returns an error if a table is
/// somehow already active (should not happen -- called at most once, at
/// startup, before config parsing).
pub fn set_active_table_from_path(path: &std::path::Path) -> Result<(), String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("failed to read unimod db `{}`: {e}", path.display()))?;
    let table = parse_csv(&text)?;
    ACTIVE_TABLE
        .set(table)
        .map_err(|_| "unimod active table already initialized".to_string())
}

fn active_table() -> &'static HashMap<u32, f32> {
    ACTIVE_TABLE.get_or_init(|| embedded_table().clone())
}

/// Resolve `UNIMOD:<id>` to its monoisotopic mass delta from the active
/// table (embedded default, or the `--unimod-db-path` override).
pub fn resolve(id: u32) -> Option<f32> {
    active_table().get(&id).copied()
}

/// Modifications resolved via `UNIMOD:<id>` *this run*, keyed by the
/// resolved mass's bit pattern -- deliberately not a full-table reverse
/// lookup. Round-tripping on output should only ever fire for a mass that
/// is bit-identical to a value this exact run resolved from the exact same
/// table, never an approximate match against an unrelated entry that
/// happens to be numerically close.
static REVERSE_TABLE: OnceLock<HashMap<u32, u32>> = OnceLock::new();

/// Set the full set of `UNIMOD:<id>` references resolved while parsing
/// this run's config, so `Peptide::Display` can print `[UNIMOD:<id>]` back
/// for any of them. One-shot (mirrors `set_active_table_from_path`) --
/// built up locally while walking the config, then set once, not mutated
/// incrementally through the global.
pub fn set_reverse_table(table: HashMap<u32, u32>) -> Result<(), String> {
    REVERSE_TABLE
        .set(table)
        .map_err(|_| "unimod reverse table already initialized".to_string())
}

/// Look up whether `mass` is bit-identical to a modification resolved via
/// `UNIMOD:<id>` this run; if so, the id to print instead of `[+mass]`.
pub fn lookup_reverse(mass: f32) -> Option<u32> {
    REVERSE_TABLE.get()?.get(&mass.to_bits()).copied()
}

/// A plain-float modification value this close to a real Unimod entry is
/// almost certainly meant to *be* that entry, imprecisely typed (e.g.
/// `57.0216` instead of Carbamidomethyl's canonical `57.021464`) --
/// treated as a config error (see `find_coincidental_match`) rather than
/// silently accepted, since it would silently fail to round-trip as
/// `[UNIMOD:<id>]` on output despite plainly being that modification.
/// Loose enough to catch a handful of rounded decimal digits, tight enough
/// that genuinely distinct modifications (which in practice differ by far
/// more than this) aren't confused for one another.
const COINCIDENCE_TOLERANCE_DA: f32 = 0.005;

/// Pure core of the coincidence check -- takes an explicit table so it's
/// testable without touching global state. Returns the closest match
/// within `tolerance`, if any (ties broken by lowest id, for determinism).
fn find_coincidental_match_in(
    table: &HashMap<u32, f32>,
    mass: f32,
    tolerance: f32,
) -> Option<(u32, f32)> {
    table
        .iter()
        .filter(|(_, &m)| (m - mass).abs() < tolerance)
        .min_by(|(id_a, m_a), (id_b, m_b)| {
            let diff_a = (*m_a - mass).abs();
            let diff_b = (*m_b - mass).abs();
            diff_a
                .partial_cmp(&diff_b)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| id_a.cmp(id_b))
        })
        .map(|(&id, &m)| (id, m))
}

/// A plain-float modification value that coincides (within
/// `COINCIDENCE_TOLERANCE_DA`) with a real entry in the active table, if
/// any -- `(unimod_id, canonical_mass)` of the closest match.
pub fn find_coincidental_match(mass: f32) -> Option<(u32, f32)> {
    find_coincidental_match_in(active_table(), mass, COINCIDENCE_TOLERANCE_DA)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn embedded_table_has_common_mods() {
        // Acetyl
        assert!((resolve(1).unwrap() - 42.010565).abs() < 1e-4);
        // Carbamidomethyl
        assert!((resolve(4).unwrap() - 57.021464).abs() < 1e-4);
    }

    #[test]
    fn unknown_id_is_none() {
        assert_eq!(resolve(u32::MAX), None);
    }

    #[test]
    fn parse_csv_rejects_malformed_line() {
        let bad = "id,name,mono_mass\nnotanumber,X,1.0\n";
        assert!(parse_csv(bad).is_err());
    }

    fn small_table() -> HashMap<u32, f32> {
        HashMap::from([(1, 42.010565), (4, 57.021464), (35, 15.994915)])
    }

    #[test]
    fn coincidental_match_catches_rounded_value() {
        // 57.0216 is what real job configs actually write for
        // Carbamidomethyl -- must be caught, not just an exact bit match.
        let found = find_coincidental_match_in(&small_table(), 57.0216, 0.005);
        assert_eq!(found, Some((4, 57.021464)));
    }

    #[test]
    fn coincidental_match_exact_value_also_caught() {
        let found = find_coincidental_match_in(&small_table(), 42.010565, 0.005);
        assert_eq!(found, Some((1, 42.010565)));
    }

    #[test]
    fn distinct_modification_not_confused() {
        // Trimethylation-ish value, far from every entry in the small
        // table -- must not match anything.
        let found = find_coincidental_match_in(&small_table(), 42.046950, 0.005);
        assert_eq!(found, None);
    }

    #[test]
    fn genuinely_novel_mass_no_match() {
        let found = find_coincidental_match_in(&small_table(), 123.456, 0.005);
        assert_eq!(found, None);
    }

    #[test]
    fn closest_match_wins_on_ties() {
        let table = HashMap::from([(1, 10.000), (2, 10.001), (3, 10.002)]);
        let found = find_coincidental_match_in(&table, 10.0015, 0.005);
        // 10.001 and 10.002 are equidistant from 10.0015 -- lower id wins.
        assert_eq!(found, Some((2, 10.001)));
    }
}
