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

    #[error("model serialisation error: {0}")]
    ModelSerialisation(#[from] ciborium::ser::Error<std::io::Error>),
}
