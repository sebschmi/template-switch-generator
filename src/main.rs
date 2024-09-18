use std::{fs::File, path::PathBuf};

use crate::error::Result;
use clap::{Args, Parser, Subcommand, ValueEnum};
use compact_genome::{
    implementation::{
        alphabets::dna_alphabet::DnaAlphabet,
        bit_array_kmer::{BitStore, BitViewSized},
        DefaultSequenceStore,
    },
    interface::{alphabet::Alphabet, sequence::GenomeSequence, sequence_store::SequenceStore},
    io::fasta::read_fasta_file,
};
use error::Error;
use n_gram_model::NGramModel;
use serde::Serialize;

mod error;
mod n_gram_model;

#[derive(Parser)]
#[command(version)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: CliCommands,
}

#[derive(Subcommand)]
enum CliCommands {
    CreateNGramModel(CreateModelCommand),
    GenerateTwin(GenerateTwinCommand),
}

#[derive(Args)]
struct CreateModelCommand {
    /// The input fasta file that contains the sequences used to create the model.
    #[arg(short, long)]
    input_fasta: PathBuf,

    /// The alphabet expected in the input file.
    #[arg(short, long, default_value = "dna")]
    alphabet: CliAlphabet,

    /// The output file in which the model is stored.
    #[arg(short, long)]
    output: PathBuf,

    /// Compute n-grams after removing all unknown characters.
    ///
    /// If this is not set, then unknown characters will result in an error.
    #[arg(short, long)]
    skip_unknown_characters: bool,

    /// Compute n-grams after capitalising all unknown characters.
    ///
    /// Capitalisation happens before skipping unknown characters.
    #[arg(short, long)]
    capitalise_characters: bool,

    /// The number of predecessor characters that determine the probability of the next character.
    ///
    /// Setting this to zero means that all characters are generated independently,
    /// setting it to one means that characters are generated with a probability depending on the previous one character, and so on.
    #[arg(short, long)]
    n_gram_context_length: usize,
}

#[derive(ValueEnum, Clone)]
enum CliAlphabet {
    Dna,
}

#[derive(Args)]
struct GenerateTwinCommand {}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        CliCommands::CreateNGramModel(create_model_command) => {
            create_n_gram_model(create_model_command)
        }
        #[expect(unused)]
        CliCommands::GenerateTwin(generate_twin_command) => todo!(),
    }
    .unwrap_or_else(|error| println!("Error: {error}"));
}

fn create_n_gram_model(create_model_command: CreateModelCommand) -> Result<()> {
    match create_model_command.alphabet {
        CliAlphabet::Dna => create_n_gram_model_for_alphabet::<{ DnaAlphabet::SIZE }, DnaAlphabet>(
            create_model_command,
        ),
    }
}

fn create_n_gram_model_for_alphabet<const ALPHABET_SIZE: usize, AlphabetType: 'static + Alphabet>(
    create_model_command: CreateModelCommand,
) -> Result<()>
where
    [u32; ALPHABET_SIZE]: Serialize,
{
    let mut sequence_store = DefaultSequenceStore::<AlphabetType>::new();
    let sequences = read_fasta_file(
        &create_model_command.input_fasta,
        &mut sequence_store,
        create_model_command.skip_unknown_characters,
        create_model_command.capitalise_characters,
    )?
    .into_iter()
    .map(|record| sequence_store.get(&record.sequence_handle));

    match create_model_command.n_gram_context_length {
        0 => create_n_gram_model_for_alphabet_and_n::<0, ALPHABET_SIZE, u8, AlphabetType, _>(
            create_model_command,
            sequences,
        ),
        1 => {
            assert!(AlphabetType::SIZE <= 256);
            create_n_gram_model_for_alphabet_and_n::<1, ALPHABET_SIZE, u8, AlphabetType, _>(
                create_model_command,
                sequences,
            )
        }
        2 => {
            assert!(AlphabetType::SIZE <= 16);
            create_n_gram_model_for_alphabet_and_n::<2, ALPHABET_SIZE, u8, AlphabetType, _>(
                create_model_command,
                sequences,
            )
        }
        3 => {
            assert!(AlphabetType::SIZE <= 4);
            create_n_gram_model_for_alphabet_and_n::<3, ALPHABET_SIZE, u8, AlphabetType, _>(
                create_model_command,
                sequences,
            )
        }
        4 => {
            assert!(AlphabetType::SIZE <= 4);
            create_n_gram_model_for_alphabet_and_n::<4, ALPHABET_SIZE, u8, AlphabetType, _>(
                create_model_command,
                sequences,
            )
        }
        n => Err(Error::UnsupportedN(n)),
    }
}

fn create_n_gram_model_for_alphabet_and_n<
    'item,
    const N: usize,
    const ALPHABET_SIZE: usize,
    BitArrayType: BitViewSized + BitStore + Serialize,
    AlphabetType: 'static + Alphabet,
    SubsequenceType: 'item + GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
>(
    create_model_command: CreateModelCommand,
    sequences: impl IntoIterator<Item = &'item SubsequenceType>,
) -> Result<()>
where
    [u32; ALPHABET_SIZE]: Serialize,
{
    let model = NGramModel::<N, ALPHABET_SIZE, _, BitArrayType>::from_sequences(sequences);
    ciborium::into_writer(&model, File::create(&create_model_command.output)?)?;

    Ok(())
}
