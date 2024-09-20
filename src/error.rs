use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("IO error: {0}")]
    IO(#[from] std::io::Error),

    #[error("genome IO error: {0}")]
    GenomeIO(#[from] compact_genome::io::error::IOError),

    #[error("n = {0} is not supported")]
    UnsupportedN(usize),

    #[error("alphabet {0:?} is not supported")]
    UnsupportedAlphabet(String),

    #[error("model serialisation error: {0}")]
    ModelSerialisation(#[from] ciborium::ser::Error<std::io::Error>),

    #[error("model deserialisation error: {0}")]
    ModelDeserialisation(#[from] ciborium::de::Error<std::io::Error>),

    #[error("the given reference ancestry fraction is not a number")]
    ReferenceAncestryFractionIsNaN,

    #[error("the given reference ancestry fraction {0} is out of range [0.0, 1.0]")]
    ReferenceAncestryFractionOutOfRange(f64),

    #[error("the given minimum template switch offset {min} is greater than the maximum {max}")]
    TemplateSwitchOffsetEmpty { min: isize, max: isize },

    #[error("the given minimum template switch length {min} is greater than the maximum {max}")]
    TemplateSwitchLengthEmpty { min: usize, max: usize },

    #[error("the given minimum template switch length difference {min} is greater than the maximum {max}")]
    TemplateSwitchLengthDifferenceEmpty { min: isize, max: isize },

    #[error("the given gap length mean is not a number")]
    GapLengthMeanIsNaN,

    #[error("the given gap length mean {actual} is out of range [{minimum}, {maximum}]")]
    GapLengthMeanOutOfRange {
        actual: f64,
        minimum: f64,
        maximum: f64,
    },

    #[error("the lambda {lambda} produced by dividing 1.0 / mean, where mean = {mean} is unsuitable for the exponential distribution")]
    GapLengthMeanLambda { lambda: f64, mean: f64 },

    #[error("the given ancestor length {length} is lower than n {n}")]
    LengthLowerThanN { length: usize, n: usize },

    #[error("the model is empty")]
    EmptyModel,

    #[error("during modifying a sequence, it became empty")]
    SequenceBecameEmpty,

    #[error("tried to generate a template switch that requires a sequence length of {template_switch_required_sequence_length} which is longer than the sequence length {sequence_length}")]
    SequenceTooShortForTemplateSwitch {
        sequence_length: usize,
        template_switch_required_sequence_length: usize,
    },

    #[error("tried to generate a gap of length {gap_length} which is longer than the sequence of length {sequence_length}")]
    SequenceTooShortForGap {
        sequence_length: usize,
        gap_length: usize,
    },
}
