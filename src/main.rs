use std::path::PathBuf;

use crate::error::Result;
use clap::{Args, Parser, Subcommand, ValueEnum};
use compact_genome::{
    implementation::{alphabets::dna_alphabet::DnaAlphabet, DefaultSequenceStore},
    interface::alphabet::Alphabet,
    io::fasta::read_fasta_file,
};

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

    /// Ignore n-grams that contain unknown characters.
    #[arg(short, long, default_value = "true")]
    ignore_unknown_characters: bool,

    /// The length of the n-grams in the model.
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

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        CliCommands::CreateNGramModel(create_model_command) => {
            create_n_gram_model(create_model_command)
        }
        #[expect(unused)]
        CliCommands::GenerateTwin(generate_twin_command) => todo!(),
    }
}

fn create_n_gram_model(create_model_command: CreateModelCommand) -> Result<()> {
    match create_model_command.alphabet {
        CliAlphabet::Dna => create_n_gram_model_for_alphabet::<DnaAlphabet>(create_model_command),
    }
}

fn create_n_gram_model_for_alphabet<AlphabetType: 'static + Alphabet>(
    create_model_command: CreateModelCommand,
) -> Result<()> {
    let mut sequence_store = DefaultSequenceStore::<AlphabetType>::new();
    #[expect(unused)]
    let handles = read_fasta_file(&create_model_command.input_fasta, &mut sequence_store)?;

    todo!()
}

#[expect(unused, clippy::extra_unused_type_parameters)]
fn create_n_gram_model_for_alphabet_and_n<const N: usize, AlphabetType: 'static + Alphabet>(
    create_model_command: CreateModelCommand,
) -> Result<()> {
    todo!()
}
