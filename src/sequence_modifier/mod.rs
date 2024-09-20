use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::{EditableGenomeSequence, GenomeSequence},
};
use rand::{seq::IteratorRandom, Rng};
use rand_distr::{Distribution, Exp};

use crate::{
    cli::{SequenceModificationAmount, SequenceModificationParameters},
    error::{Error, Result},
};

pub struct SequenceModifier {
    sequence_modification_amount: SequenceModificationAmount,
    sequence_modification_parameters: SequenceModificationParameters,
}

#[derive(Clone, Copy)]
pub enum SequenceModification {
    #[expect(unused)]
    TemplateSwitch {
        position: usize,
        length: usize,
        offset: isize,
        length_difference: isize,
    },
    Insertion {
        position: usize,
        source: usize,
        length: usize,
    },
    Deletion {
        position: usize,
        length: usize,
    },
    Substitution {
        position: usize,
        character_increment: usize,
    },
}

pub struct SequenceModifierPair {
    pub reference_modifier: SequenceModifier,
    pub query_modifier: SequenceModifier,
}

impl SequenceModifier {
    pub fn new_modifier_pair(
        reference_ancestry_fraction: f64,
        sequence_modification_amount: SequenceModificationAmount,
        sequence_modification_parameters: SequenceModificationParameters,
        rng: &mut impl Rng,
    ) -> SequenceModifierPair {
        let (query_template_switch_amount, reference_template_switch_amount) = split_int_random(
            sequence_modification_amount.template_switch_amount,
            reference_ancestry_fraction,
            rng,
        );
        let (query_gap_amount, reference_gap_amount) = split_int_random(
            sequence_modification_amount.gap_amount,
            reference_ancestry_fraction,
            rng,
        );
        let (query_substitution_amount, reference_substitution_amount) = split_int_random(
            sequence_modification_amount.substitution_amount,
            reference_ancestry_fraction,
            rng,
        );

        SequenceModifierPair {
            reference_modifier: SequenceModifier {
                sequence_modification_amount: SequenceModificationAmount {
                    template_switch_amount: reference_template_switch_amount,
                    gap_amount: reference_gap_amount,
                    substitution_amount: reference_substitution_amount,
                },
                sequence_modification_parameters,
            },
            query_modifier: SequenceModifier {
                sequence_modification_amount: SequenceModificationAmount {
                    template_switch_amount: query_template_switch_amount,
                    gap_amount: query_gap_amount,
                    substitution_amount: query_substitution_amount,
                },
                sequence_modification_parameters,
            },
        }
    }

    pub fn next(
        &mut self,
        sequence_length: usize,
        alphabet_size: usize,
        rng: &mut impl Rng,
    ) -> Result<Option<SequenceModification>> {
        assert!(alphabet_size > 1);
        if sequence_length == 0 {
            return Err(Error::SequenceBecameEmpty);
        }

        let sum = self.sequence_modification_amount.template_switch_amount
            + self.sequence_modification_amount.gap_amount
            + self.sequence_modification_amount.substitution_amount;

        if let Some(index) = (0..sum).choose(rng) {
            Ok(Some(
                if index < self.sequence_modification_amount.template_switch_amount {
                    debug_assert!(self.sequence_modification_amount.template_switch_amount > 0);
                    self.sequence_modification_amount.template_switch_amount -= 1;

                    let offset = (self
                        .sequence_modification_parameters
                        .template_switch_min_offset
                        ..=self
                            .sequence_modification_parameters
                            .template_switch_max_offset)
                        .choose(rng)
                        .unwrap();
                    let length = (self
                        .sequence_modification_parameters
                        .template_switch_min_length
                        ..=self
                            .sequence_modification_parameters
                            .template_switch_max_length)
                        .choose(rng)
                        .unwrap() as isize;
                    let length_difference = (self
                        .sequence_modification_parameters
                        .template_switch_min_length_difference
                        ..=self
                            .sequence_modification_parameters
                            .template_switch_max_length_difference)
                        .choose(rng)
                        .unwrap();
                    let position_range = 0.max(offset - length)
                        + self.sequence_modification_parameters.template_switch_margin as isize
                        ..(sequence_length as isize
                            - 0.max(offset).max(length).max(length + length_difference))
                            - self.sequence_modification_parameters.template_switch_margin as isize;
                    #[expect(unused)]
                    let position = position_range.clone().choose(rng).ok_or_else(|| {
                        Error::SequenceTooShortForTemplateSwitch {
                            sequence_length,
                            template_switch_required_sequence_length: (position_range.start
                                + (sequence_length as isize - position_range.end))
                                as usize,
                        }
                    })?;

                    /*SequenceModification::TemplateSwitch {
                        position: position as usize,
                        length: length as usize,
                        offset,
                        length_difference,
                    }*/
                    todo!("prevent consecutive template switches from interacting")
                } else if index
                    < self.sequence_modification_amount.template_switch_amount
                        + self.sequence_modification_amount.gap_amount
                {
                    debug_assert!(self.sequence_modification_amount.gap_amount > 0);
                    self.sequence_modification_amount.gap_amount -= 1;

                    let mean_gap_length = self.sequence_modification_parameters.gap_length_mean;
                    let gap_length = Exp::new(1.0 / mean_gap_length)
                        .map_err(|error| match error {
                            rand_distr::ExpError::LambdaTooSmall => Error::GapLengthMeanLambda {
                                lambda: 1.0 / mean_gap_length,
                                mean: mean_gap_length,
                            },
                        })?
                        .sample(rng);
                    let gap_length = if gap_length < 1.0 {
                        1
                    } else {
                        gap_length.round() as usize
                    };

                    if gap_length < sequence_length {
                        return Err(Error::SequenceTooShortForGap {
                            sequence_length,
                            gap_length,
                        });
                    }

                    if rng.gen_bool(0.5) {
                        SequenceModification::Insertion {
                            position: (0..sequence_length).choose(rng).unwrap(),
                            source: (0..sequence_length - gap_length).choose(rng).unwrap(),
                            length: gap_length,
                        }
                    } else {
                        SequenceModification::Deletion {
                            position: (0..sequence_length - gap_length).choose(rng).unwrap(),
                            length: gap_length,
                        }
                    }
                } else {
                    debug_assert!(self.sequence_modification_amount.substitution_amount > 0);
                    self.sequence_modification_amount.substitution_amount -= 1;

                    SequenceModification::Substitution {
                        position: (0..sequence_length).choose(rng).unwrap(),
                        character_increment: (1..alphabet_size).choose(rng).unwrap(),
                    }
                },
            ))
        } else {
            Ok(None)
        }
    }

    pub fn apply<
        AlphabetType: Alphabet,
        SequenceType: EditableGenomeSequence<AlphabetType, SubsequenceType>,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &mut self,
        sequence: &mut SequenceType,
        rng: &mut impl Rng,
    ) -> Result<()> {
        while let Some(modification) = self.next(sequence.len(), AlphabetType::SIZE, rng)? {
            modification.apply(sequence)?;
        }

        Ok(())
    }
}

impl SequenceModification {
    #[expect(unused)]
    pub fn apply<
        AlphabetType: Alphabet,
        SequenceType: EditableGenomeSequence<AlphabetType, SubsequenceType>,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        sequence: &mut SequenceType,
    ) -> Result<()> {
        match *self {
            SequenceModification::TemplateSwitch {
                position,
                length,
                offset,
                length_difference,
            } => todo!(),
            SequenceModification::Insertion {
                position,
                source,
                length,
            } => todo!(),
            SequenceModification::Deletion { position, length } => todo!(),
            SequenceModification::Substitution {
                position,
                character_increment,
            } => {
                assert!(character_increment > 0);
                let character_index = sequence[position].index();
                let character_index = (character_index + character_increment) % AlphabetType::SIZE;
                sequence.set(
                    position,
                    AlphabetType::CharacterType::from_index(character_index).unwrap(),
                );
            }
        }

        Ok(())
    }
}

fn split_int_random(int: usize, fraction: f64, rng: &mut impl Rng) -> (usize, usize) {
    assert!(fraction >= 0.0);
    assert!(fraction <= 1.0);
    assert!(fraction.is_finite());

    let amount1 = fraction * int as f64;
    let amount2 = int as f64 - amount1;
    let amount1_int = amount1.floor() as usize;
    let amount2_int = amount2.floor() as usize;

    if amount1_int + amount2_int == int - 1 {
        if rng.gen_bool(fraction) {
            (amount1_int + 1, amount2_int)
        } else {
            (amount1_int, amount2_int + 1)
        }
    } else {
        assert_eq!(amount1_int + amount2_int, int);
        (amount1_int, amount2_int)
    }
}
