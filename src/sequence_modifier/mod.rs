use compact_genome::interface::{
    alphabet::{Alphabet, AlphabetCharacter},
    sequence::{EditableGenomeSequence, GenomeSequence},
};
use rand::{seq::IteratorRandom, Rng};
use rand_distr::{Distribution, Exp};
use template_switch_overlap_detector::TemplateSwitchOverlapDetector;

use crate::{
    cli::{SequenceModificationAmount, SequenceModificationParameters},
    error::{Error, Result},
};

pub mod template_switch_overlap_detector;

pub struct SequenceModifier {
    sequence_modification_amount: SequenceModificationAmount,
    sequence_modification_parameters: SequenceModificationParameters,
}

#[derive(Debug, Clone, Copy)]
pub enum SequenceModification {
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
        template_switch_overlap_detector: &mut TemplateSwitchOverlapDetector,
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
                    let mut tries = 0;

                    loop {
                        if tries
                            < self
                                .sequence_modification_parameters
                                .template_switch_maximum_overlap_tries
                        {
                            tries += 1;
                        } else {
                            return Err(Error::TemplateSwitchOverlap);
                        }

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
                            ..=(self
                                .sequence_modification_parameters
                                .template_switch_max_length_difference)
                                .min(length))
                            .choose(rng)
                            .unwrap();
                        let position_range = 0.max(offset - length)
                            + self.sequence_modification_parameters.template_switch_margin as isize
                            ..(sequence_length as isize
                                - 0.max(offset).max(length).max(length + length_difference))
                                - self.sequence_modification_parameters.template_switch_margin
                                    as isize;

                        let position = position_range.clone().choose(rng).ok_or_else(|| {
                            Error::SequenceTooShortForTemplateSwitch {
                                sequence_length,
                                template_switch_required_sequence_length: (position_range.start
                                    + (sequence_length as isize - position_range.end))
                                    as usize,
                            }
                        })?;

                        let result = SequenceModification::TemplateSwitch {
                            position: position as usize,
                            length: length as usize,
                            offset,
                            length_difference,
                        };

                        if self
                            .sequence_modification_parameters
                            .template_switch_overlap
                        {
                            break result;
                        } else {
                            match template_switch_overlap_detector
                                .apply_modification(result) {
                                    template_switch_overlap_detector::TemplateSwitchCollision::Overlap => { /* retry */ }
                                    template_switch_overlap_detector::TemplateSwitchCollision::Independent => break result,
                                }
                        }
                    }
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

                    if gap_length > sequence_length {
                        return Err(Error::SequenceTooShortForGap {
                            sequence_length,
                            gap_length,
                        });
                    }

                    let result = if rng.gen_bool(0.5) {
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
                    };

                    if !self
                        .sequence_modification_parameters
                        .template_switch_overlap
                    {
                        template_switch_overlap_detector.apply_modification(result);
                    }
                    result
                } else {
                    debug_assert!(self.sequence_modification_amount.substitution_amount > 0);
                    self.sequence_modification_amount.substitution_amount -= 1;

                    let result = SequenceModification::Substitution {
                        position: (0..sequence_length).choose(rng).unwrap(),
                        character_increment: (1..alphabet_size).choose(rng).unwrap(),
                    };

                    if !self
                        .sequence_modification_parameters
                        .template_switch_overlap
                    {
                        template_switch_overlap_detector.apply_modification(result);
                    }
                    result
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
        template_switch_overlap_detector: &mut TemplateSwitchOverlapDetector,
        rng: &mut impl Rng,
    ) -> Result<()> {
        while let Some(modification) = self.next(
            sequence.len(),
            AlphabetType::SIZE,
            template_switch_overlap_detector,
            rng,
        )? {
            modification.apply(sequence)?;
        }

        Ok(())
    }
}

impl SequenceModification {
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
            } => {
                let replacement: SequenceType = sequence
                    .reverse_complement_iter()
                    .skip(sequence.len() - (position as isize + offset + 1) as usize)
                    .take(length)
                    .collect();
                sequence.splice(
                    position..((position as isize + length as isize - length_difference) as usize),
                    replacement,
                );
            }

            SequenceModification::Insertion {
                position,
                source,
                length,
            } => {
                let insertion: SequenceType =
                    sequence.iter().skip(source).take(length).cloned().collect();
                sequence.splice(position..position, insertion);
            }

            SequenceModification::Deletion { position, length } => {
                sequence.delete(position..position + length)
            }

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
