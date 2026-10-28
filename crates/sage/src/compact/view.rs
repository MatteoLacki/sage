//! A [`CompactPeptides`] row read as a [`PeptideView`], without building a
//! [`Peptide`](crate::peptide::Peptide): residues straight from the proteome copy,
//! modifications decided exactly as [`CompactPeptides::materialize`] decides them
//! (sage's own `Peptide::apply_site`/`Peptide::static_mods`, then
//! `Peptide::reverse` for decoys). See `docs/ai/peptide_views.md`.

use super::{CompactPeptide, CompactPeptides};
use crate::enzyme::Position;
use crate::modification::ModificationSpecificity;
use crate::peptide::Site;
use crate::peptide_view::PeptideView;
use std::collections::HashMap;

/// The run's static mods in the order `materialize` applies them (its `HashMap`'s
/// iteration order), and each residue's residue rule.
#[derive(Clone, Debug)]
pub(super) struct StaticMods {
    ordered: Vec<(ModificationSpecificity, f32)>,
    residue: [Option<f32>; 256],
}

impl StaticMods {
    pub(super) fn new(static_mods: &HashMap<ModificationSpecificity, f32>) -> Self {
        let ordered: Vec<(ModificationSpecificity, f32)> = static_mods
            .iter()
            .map(|(&target, &mass)| (target, mass))
            .collect();
        let mut residue = [None; 256];
        for &(target, mass) in &ordered {
            if let ModificationSpecificity::Residue(r) = target {
                residue[r as usize] = Some(mass);
            }
        }
        Self { ordered, residue }
    }
}

/// One row of a [`CompactPeptides`] table as a [`PeptideView`].
#[derive(Clone, Copy)]
pub struct CompactView<'t> {
    table: &'t CompactPeptides,
    row: CompactPeptide,
    /// The target sequence's residues (a decoy reads them interior-reversed).
    target: &'t [u8],
}

impl CompactPeptides {
    /// Row `row` as a [`PeptideView`]: what `materialize` would build, read in place.
    pub fn view(&self, row: usize) -> CompactView<'_> {
        self.view_of(self.rows.get(row))
    }

    pub fn view_of(&self, row: CompactPeptide) -> CompactView<'_> {
        CompactView {
            table: self,
            row,
            target: self.proteome.window(row.start, row.length),
        }
    }
}

impl<'t> CompactView<'t> {
    /// The target this row reverses, read from the same row; a target itself.
    pub fn target(&self) -> CompactView<'t> {
        CompactView {
            row: self.row.target_of(),
            ..*self
        }
    }

    /// The target's residues: the proteome window the row reads.
    pub fn target_residues(&self) -> &'t [u8] {
        self.target
    }

    pub fn semi_enzymatic(&self) -> bool {
        self.row.is_semi_enzymatic()
    }

    fn last(&self) -> usize {
        self.target.len() - 1
    }

    /// Where residue `index` sits in the target: a decoy is its target with the
    /// interior reversed, the first and last residues kept (`Peptide::reverse`).
    fn target_index(&self, index: usize) -> usize {
        let last = self.last();
        match self.row.is_decoy() && last > 1 && index > 0 && index < last {
            true => last - index,
            false => index,
        }
    }

    fn variable_mass(&self, site: Site) -> Option<f32> {
        self.row.var_sites.iter().find_map(|var_site| {
            (var_site.site()? == site)
                .then(|| self.table.variable_mod_masses[var_site.mass_id as usize])
        })
    }

    fn at_protein_n(&self) -> bool {
        matches!(self.row.position(), Position::Nterm | Position::Full)
    }

    fn at_protein_c(&self) -> bool {
        matches!(self.row.position(), Position::Cterm | Position::Full)
    }

    /// Whether static rule `target` puts its mass on residue `index` of the target.
    fn static_rule_hits(&self, target: ModificationSpecificity, index: usize) -> bool {
        let last = self.last();
        match target {
            ModificationSpecificity::Residue(r) => self.target[index] == r,
            ModificationSpecificity::PeptideN(Some(r)) => index == 0 && self.target[0] == r,
            ModificationSpecificity::PeptideC(Some(r)) => index == last && self.target[last] == r,
            ModificationSpecificity::ProteinN(Some(r)) => {
                index == 0 && self.at_protein_n() && self.target[0] == r
            }
            ModificationSpecificity::ProteinC(Some(r)) => {
                index == last && self.at_protein_c() && self.target[last] == r
            }
            ModificationSpecificity::PeptideN(None)
            | ModificationSpecificity::PeptideC(None)
            | ModificationSpecificity::ProteinN(None)
            | ModificationSpecificity::ProteinC(None) => false,
        }
    }

    /// The mass `materialize` leaves on residue `index` of the target. A residue takes
    /// a mass only while it still holds 0: a variable site there first, then the
    /// static rules in their order (`apply_site` adds to the 0, a residue rule
    /// assigns).
    fn target_modification(&self, index: usize) -> f32 {
        if let Some(mass) = self.variable_mass(Site::Sequence(index as u32)) {
            let applied = 0.0 + mass;
            if applied != 0.0 {
                return applied;
            }
        }
        let interior = index != 0 && index != self.last();
        if interior {
            return self.table.static_mods.residue[self.target[index] as usize].unwrap_or(0.0);
        }
        for &(target, mass) in &self.table.static_mods.ordered {
            if !self.static_rule_hits(target, index) {
                continue;
            }
            let applied = match target {
                ModificationSpecificity::Residue(_) => mass,
                _ => 0.0 + mass,
            };
            if applied != 0.0 {
                return applied;
            }
        }
        0.0
    }

    /// A terminus' modification: the variable site there, else the first static rule
    /// for it. `apply_site` sets a terminus only while it is `None`.
    fn terminal(
        &self,
        site: Site,
        takes_rule: impl Fn(ModificationSpecificity) -> bool,
    ) -> Option<f32> {
        if let Some(mass) = self.variable_mass(site) {
            return Some(0.0 + mass);
        }
        self.table
            .static_mods
            .ordered
            .iter()
            .find(|&&(target, _)| takes_rule(target))
            .map(|&(_, mass)| 0.0 + mass)
    }
}

impl PeptideView for CompactView<'_> {
    fn len(&self) -> usize {
        self.target.len()
    }

    fn residue(&self, index: usize) -> u8 {
        self.target[self.target_index(index)]
    }

    fn modification(&self, index: usize) -> f32 {
        self.target_modification(self.target_index(index))
    }

    fn nterm(&self) -> Option<f32> {
        self.terminal(Site::Nterm, |target| match target {
            ModificationSpecificity::PeptideN(None) => true,
            ModificationSpecificity::ProteinN(None) => self.at_protein_n(),
            _ => false,
        })
    }

    fn cterm(&self) -> Option<f32> {
        self.terminal(Site::Cterm, |target| match target {
            ModificationSpecificity::PeptideC(None) => true,
            ModificationSpecificity::ProteinC(None) => self.at_protein_c(),
            _ => false,
        })
    }

    fn monoisotopic(&self) -> f32 {
        self.row.monoisotopic
    }

    fn decoy(&self) -> bool {
        self.row.is_decoy()
    }

    fn missed_cleavages(&self) -> u8 {
        self.row.missed_cleavages()
    }
}
