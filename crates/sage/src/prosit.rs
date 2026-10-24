//! Prosit 2023 timsTOF inputs, read directly from a peptide view.
//! The vocabulary is specific to this model, not a general modification encoder.
use crate::peptide_view::PeptideView;
use std::collections::HashMap;

pub const SEQUENCE_LENGTH: usize = 30;
pub const FRAGMENTS: usize = 174;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModificationPolicy {
    Original,
    StripAll,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Unsupported {
    Length,
    Residue,
    TerminalModification,
    ResidueModification,
}

impl Unsupported {
    pub fn reason(self) -> &'static str {
        match self {
            Self::Length => "length",
            Self::Residue => "residue",
            Self::TerminalModification => "terminal_modification",
            Self::ResidueModification => "residue_modification",
        }
    }
}

pub struct Tokenizer {
    policy: ModificationPolicy,
    // Only modifications identified by the run's config are admitted. No fuzzy
    // mass matching: different modifications may have the same mass.
    modifications: HashMap<u32, u32>,
}

impl Tokenizer {
    pub fn new(policy: ModificationPolicy, modifications: HashMap<u32, u32>) -> Self {
        Self {
            policy,
            modifications,
        }
    }

    pub fn encode(
        &self,
        peptide: &impl PeptideView,
        output: &mut [i32; SEQUENCE_LENGTH],
    ) -> Result<(), Unsupported> {
        output.fill(0);
        if peptide.is_empty() || peptide.len() > SEQUENCE_LENGTH {
            return Err(Unsupported::Length);
        }
        if self.policy == ModificationPolicy::Original
            && (peptide.nterm().is_some() || peptide.cterm().is_some())
        {
            return Err(Unsupported::TerminalModification);
        }
        for (index, token) in output.iter_mut().enumerate().take(peptide.len()) {
            let residue = peptide.residue(index);
            *token = match residue {
                b'A' => 1,
                b'C' => 2,
                b'D' => 3,
                b'E' => 4,
                b'F' => 5,
                b'G' => 6,
                b'H' => 7,
                b'I' => 8,
                b'K' => 9,
                b'L' => 10,
                b'M' => 11,
                b'N' => 12,
                b'P' => 13,
                b'Q' => 14,
                b'R' => 15,
                b'S' => 16,
                b'T' => 17,
                b'V' => 18,
                b'W' => 19,
                b'Y' => 20,
                _ => return Err(Unsupported::Residue),
            };
            if self.policy == ModificationPolicy::StripAll {
                continue;
            }
            let mass = peptide.modification(index);
            if mass == 0.0 {
                continue;
            }
            match (residue, self.modifications.get(&mass.to_bits())) {
                (b'M', Some(35)) => *token = 21,
                (b'C', Some(4)) | (b'K', Some(259)) | (b'R', Some(267)) => (),
                _ => return Err(Unsupported::ResidueModification),
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{compact::CompactPeptides, database::Builder, fasta::Fasta, peptide::Peptide};

    fn tokenizer(policy: ModificationPolicy) -> Tokenizer {
        let ids = [4, 35, 259, 267];
        Tokenizer::new(
            policy,
            ids.into_iter()
                .map(|id| (crate::unimod::resolve(id).unwrap().to_bits(), id))
                .collect(),
        )
    }

    fn peptide(sequence: &[u8]) -> Peptide {
        // Test through the real digestion constructor, including modification arrays.
        let fasta = Fasta::parse(
            format!(">p\n{}\n", std::str::from_utf8(sequence).unwrap()),
            "rev_",
            true,
        );
        let mut builder: Builder = serde_json::from_value(serde_json::json!({
            "enzyme": {"cleave_at": "", "restrict": null, "min_len": sequence.len(), "max_len": sequence.len()},
            "peptide_min_mass": 0.0, "peptide_max_mass": 100000.0,
            "static_mods": {}, "variable_mods": {}
        })).unwrap();
        builder.fasta = Some("fixture".into());
        builder
            .make_parameters()
            .digest(&fasta)
            .into_iter()
            .find(|p| !p.decoy)
            .unwrap()
    }

    #[test]
    fn exact_vocabulary_padding_and_modified_tokens() {
        let mut p = peptide(b"ACDMKR");
        for (position, id) in [(1, 4), (3, 35), (4, 259), (5, 267)] {
            p.modifications[position] = crate::unimod::resolve(id).unwrap();
        }
        let mut tokens = [99; SEQUENCE_LENGTH];
        tokenizer(ModificationPolicy::Original)
            .encode(&p, &mut tokens)
            .unwrap();
        assert_eq!(&tokens[..6], &[1, 2, 3, 21, 9, 15]);
        assert!(tokens[6..].iter().all(|&t| t == 0));
        tokenizer(ModificationPolicy::StripAll)
            .encode(&p, &mut tokens)
            .unwrap();
        assert_eq!(&tokens[..6], &[1, 2, 3, 11, 9, 15]);
    }

    #[test]
    fn rejects_unrepresented_modifications_and_lengths() {
        let mut p = peptide(b"ACDMK");
        let mut tokens = [0; SEQUENCE_LENGTH];
        p.modifications[2] = 15.994915;
        assert_eq!(
            tokenizer(ModificationPolicy::Original).encode(&p, &mut tokens),
            Err(Unsupported::ResidueModification)
        );
        p.nterm = Some(42.010565);
        assert_eq!(
            tokenizer(ModificationPolicy::Original).encode(&p, &mut tokens),
            Err(Unsupported::TerminalModification)
        );
        tokenizer(ModificationPolicy::StripAll)
            .encode(&p, &mut tokens)
            .unwrap();
        p.sequence = std::sync::Arc::from(&b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"[..]);
        assert_eq!(
            tokenizer(ModificationPolicy::StripAll).encode(&p, &mut tokens),
            Err(Unsupported::Length)
        );
    }

    #[test]
    fn compact_and_materialized_tokens_match_with_decoys_and_terminal_rules() {
        let fasta = Fasta::parse(">p\nMACDMKRACDMK\n".into(), "rev_", true);
        let mut builder: Builder = serde_json::from_value(serde_json::json!({
            "enzyme": {"cleave_at": "", "restrict": null, "min_len": 5, "max_len": 9},
            "peptide_min_mass": 0.0, "peptide_max_mass": 100000.0,
            "static_mods": {"C": crate::unimod::resolve(4).unwrap()},
            "variable_mods": {"M": [crate::unimod::resolve(35).unwrap()], "[": [42.010565]},
            "max_variable_mods": 2
        }))
        .unwrap();
        builder.fasta = Some("fixture".into());
        let table = CompactPeptides::build(&builder.make_parameters(), &fasta).unwrap();
        assert!(table.rows.len() > 0);
        let mut decoys = 0;
        for row in 0..table.rows.len() {
            let view = table.view(row);
            if view.decoy() {
                decoys += 1;
            }
            let p = table.materialize(&table.rows.get(row));
            for policy in [ModificationPolicy::Original, ModificationPolicy::StripAll] {
                let tokenizer = tokenizer(policy);
                let (mut a, mut b) = ([0; SEQUENCE_LENGTH], [0; SEQUENCE_LENGTH]);
                assert_eq!(
                    tokenizer.encode(&view, &mut a),
                    tokenizer.encode(&p, &mut b)
                );
                assert_eq!(a, b);
            }
        }
        assert!(decoys > 0);
    }
}
