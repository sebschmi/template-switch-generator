use thiserror::Error;

pub type Result<T> = std::result::Result<T, Error>;

#[derive(Error, Debug)]
pub enum Error {
    #[error("Genome IO error: {0}")]
    GenomeIOError(#[from] compact_genome::io::error::IOError),
}
