use std::{io::Write, ops::Range};

use crate::{cli::SequenceModificationParameters, error::Result};

use super::SequenceModification;

#[derive(Debug)]
pub struct TemplateSwitchOverlapDetector {
    template_switches: Vec<Range<usize>>,
    modification_stack: Vec<SequenceModification>,
    margin: usize,
}

#[derive(Debug, Eq, PartialEq)]
pub enum TemplateSwitchCollision {
    Overlap,
    Independent,
}

impl TemplateSwitchOverlapDetector {
    pub fn new(sequence_modification_parameters: &SequenceModificationParameters) -> Self {
        Self::from_template_switch_margin(sequence_modification_parameters.template_switch_margin)
    }

    fn from_template_switch_margin(margin: usize) -> Self {
        Self {
            template_switches: Default::default(),
            modification_stack: Default::default(),
            margin,
        }
    }

    pub fn clear_modification_stack(&mut self) {
        self.modification_stack.clear();
    }

    pub fn apply_modification(
        &mut self,
        sequence_modification: SequenceModification,
    ) -> TemplateSwitchCollision {
        match sequence_modification {
            SequenceModification::TemplateSwitch {
                position,
                length,
                offset,
                length_difference,
            } => {
                let range_offset =
                    position.min((position as isize - length as isize + offset) as usize);
                let range_limit = position
                    .max((position as isize + length as isize - length_difference) as usize)
                    .max((position as isize + offset) as usize);

                debug_assert!(range_offset <= range_limit);
                if range_offset < self.margin
                    || range_offset > isize::MAX as usize
                    || range_limit > isize::MAX as usize
                {
                    return TemplateSwitchCollision::Overlap;
                }

                let new_range = (range_offset - self.margin)..(range_limit + self.margin);

                let new_range = self.modification_stack.iter().rev().fold(
                    new_range,
                    |mut new_range, sequence_modification| match *sequence_modification {
                        SequenceModification::TemplateSwitch {
                            position,
                            length_difference,
                            ..
                        } => {
                            if new_range.start > position {
                                new_range.start = position
                                    .max((new_range.start as isize - length_difference) as usize);
                            }
                            if new_range.end > position {
                                new_range.end = position
                                    .max((new_range.end as isize - length_difference) as usize);
                            }
                            new_range
                        }

                        SequenceModification::Insertion {
                            position, length, ..
                        } => {
                            if new_range.start > position {
                                new_range.start = position.max(new_range.start - length);
                            }
                            if new_range.end > position {
                                new_range.end = position.max(new_range.end - length);
                            }
                            new_range
                        }

                        SequenceModification::Deletion { position, length } => {
                            if new_range.start > position {
                                new_range.start = position.max(new_range.start + length);
                            }
                            if new_range.end > position {
                                new_range.end = position.max(new_range.end + length);
                            }
                            new_range
                        }

                        SequenceModification::Substitution { .. } => {
                            // No range modification.
                            new_range
                        }
                    },
                );

                let insertion_offset = self
                    .template_switches
                    .iter()
                    .take_while(|range| range.end <= new_range.start)
                    .count();
                if let Some(range) = self.template_switches.get(insertion_offset) {
                    if new_range.start < range.end && range.start < new_range.end {
                        TemplateSwitchCollision::Overlap
                    } else {
                        self.template_switches
                            .insert(insertion_offset, new_range.clone());
                        self.modification_stack.push(sequence_modification);
                        TemplateSwitchCollision::Independent
                    }
                } else {
                    self.template_switches
                        .insert(insertion_offset, new_range.clone());
                    self.modification_stack.push(sequence_modification);
                    TemplateSwitchCollision::Independent
                }
            }

            SequenceModification::Insertion { .. }
            | SequenceModification::Deletion { .. }
            | SequenceModification::Substitution { .. } => {
                self.modification_stack.push(sequence_modification);

                TemplateSwitchCollision::Independent
            }
        }
    }

    pub fn write_modifications(&self, output: &mut impl Write) -> Result<()> {
        for modification in &self.modification_stack {
            writeln!(output, "{modification}")?;
        }

        Ok(())
    }
}

#[cfg(test)]
#[allow(clippy::single_range_in_vec_init)]
mod tests {
    use crate::sequence_modifier::{
        template_switch_overlap_detector::TemplateSwitchCollision, SequenceModification,
    };

    use super::TemplateSwitchOverlapDetector;

    #[test]
    fn simple() {
        let mut tsod = TemplateSwitchOverlapDetector::from_template_switch_margin(10);
        assert_eq!(
            tsod.apply_modification(SequenceModification::TemplateSwitch {
                position: 50,
                length: 10,
                offset: -5,
                length_difference: 5,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(tsod.template_switches.as_slice(), [25..65]);
        assert_eq!(
            tsod.apply_modification(SequenceModification::TemplateSwitch {
                position: 100,
                length: 10,
                offset: -5,
                length_difference: 5,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(tsod.template_switches.as_slice(), [25..65, 70..110]);
        assert_eq!(
            tsod.apply_modification(SequenceModification::TemplateSwitch {
                position: 150,
                length: 10,
                offset: -5,
                length_difference: -10,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.template_switches.as_slice(),
            [25..65, 70..110, 115..170]
        );

        assert_eq!(
            tsod.apply_modification(SequenceModification::Substitution {
                position: 60,
                character_increment: 2
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.apply_modification(SequenceModification::Substitution {
                position: 70,
                character_increment: 2
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.apply_modification(SequenceModification::Substitution {
                position: 80,
                character_increment: 2
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.apply_modification(SequenceModification::TemplateSwitch {
                position: 200,
                length: 10,
                offset: 20,
                length_difference: -1,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.template_switches.as_slice(),
            [25..65, 70..110, 115..170, 190..230]
        );

        assert_eq!(
            tsod.apply_modification(SequenceModification::Insertion {
                position: 70,
                source: 10,
                length: 1,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.apply_modification(SequenceModification::Deletion {
                position: 70,
                length: 10,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.apply_modification(SequenceModification::TemplateSwitch {
                position: 250,
                length: 20,
                offset: 10,
                length_difference: 0,
            }),
            TemplateSwitchCollision::Independent
        );
        assert_eq!(
            tsod.template_switches.as_slice(),
            [25..65, 70..110, 115..170, 190..230, 240..290]
        );

        assert_eq!(
            tsod.apply_modification(SequenceModification::TemplateSwitch {
                position: 200,
                length: 20,
                offset: 10,
                length_difference: 0,
            }),
            TemplateSwitchCollision::Overlap
        );
        assert_eq!(
            tsod.template_switches.as_slice(),
            [25..65, 70..110, 115..170, 190..230, 240..290]
        );
    }
}
