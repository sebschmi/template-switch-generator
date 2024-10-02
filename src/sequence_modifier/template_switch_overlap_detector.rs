use std::ops::Range;

use crate::cli::SequenceModificationParameters;

use super::SequenceModification;

#[derive(Debug)]
pub struct TemplateSwitchOverlapDetector {
    template_switches: Vec<Range<usize>>,
    modification_stack: Vec<SequenceModification>,
    margin: usize,
}

pub enum TemplateSwitchCollision {
    Overlap,
    Independent,
}

impl TemplateSwitchOverlapDetector {
    pub fn new(sequence_modification_parameters: &SequenceModificationParameters) -> Self {
        Self {
            template_switches: Default::default(),
            modification_stack: Default::default(),
            margin: sequence_modification_parameters.template_switch_margin,
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
                    .max((position as isize + length_difference) as usize)
                    .max((position as isize + offset) as usize);
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
                    .skip_while(|range| range.end <= new_range.start)
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
}
