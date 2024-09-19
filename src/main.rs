use std::{
    fs::File,
    io::{BufReader, BufWriter},
};

use crate::error::Result;
use choose_alphabet_and_n::{call, ChooseAlphabetAndN};
use clap::{Parser, ValueEnum};
use cli::{
    Cli, CliAlphabet, CliCommands, CreateModelCommand, GeneratePairCommand, IntoCliAlphabet,
};
use compact_genome::{
    implementation::{
        bit_array_kmer::{BitStore, BitView, BitViewSized},
        handle_sequence_store::HandleSequenceStore,
        DefaultGenome, DefaultSubGenome,
    },
    interface::alphabet::Alphabet,
    io::fasta::{read_fasta_file, write_fasta_file, FastaRecord},
};
use error::Error;
use n_gram_model::NGramModel;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use serde::{Deserialize, Serialize};

mod choose_alphabet_and_n;
mod cli;
mod error;
mod n_gram_model;

fn main() {
    let cli = Cli::parse();

    match cli.command {
        CliCommands::CreateNGramModel(create_model_command) => call::<CreateNGramModel>(
            create_model_command.alphabet,
            create_model_command.n_gram_context_length,
            create_model_command,
        ),
        CliCommands::GeneratePair(generate_pair_command) => generate_pair(generate_pair_command),
    }
    .unwrap_or_else(|error| println!("Error: {error}"));
}

struct CreateNGramModel;

impl ChooseAlphabetAndN for CreateNGramModel {
    type Arguments = CreateModelCommand;

    type Return = ();

    fn call<
        const N: usize,
        const ALPHABET_SIZE: usize,
        BitArrayType: BitViewSized + BitStore + Serialize + for<'de> Deserialize<'de>,
        AlphabetType: 'static + Alphabet + IntoCliAlphabet,
    >(
        create_model_command: Self::Arguments,
    ) -> Result<Self::Return>
    where
        [u32; ALPHABET_SIZE]: Serialize + for<'de> Deserialize<'de>,
    {
        // Load sequences.
        let mut sequence_store =
            HandleSequenceStore::<AlphabetType, DefaultGenome<_>, DefaultSubGenome<_>>::new();
        let sequences = read_fasta_file(
            &create_model_command.input_fasta,
            &mut sequence_store,
            create_model_command.skip_unknown_characters,
            create_model_command.capitalise_characters,
        )?
        .into_iter()
        .map(|record| record.sequence_handle);

        // Create model.
        let model = NGramModel::<N, ALPHABET_SIZE, _, BitArrayType>::from_sequences(sequences);
        let mut output = BufWriter::new(File::create(&create_model_command.output)?);

        // Write model parameters and model.
        ciborium::into_writer(&N, &mut output)?;
        ciborium::into_writer(
            &AlphabetType::into_cli_alphabet()
                .to_possible_value()
                .unwrap()
                .get_name(),
            &mut output,
        )?;
        ciborium::into_writer(&model, &mut output)?;

        Ok(())
    }
}

fn generate_pair(generate_pair_command: GeneratePairCommand) -> Result<()> {
    // Verify some inputs.
    if generate_pair_command.reference_ancestry_fraction.is_nan() {
        return Err(Error::ReferenceAncestryFractionIsNaN);
    }
    if generate_pair_command.reference_ancestry_fraction < 0.0
        || generate_pair_command.reference_ancestry_fraction > 1.0
    {
        return Err(Error::ReferenceAncestryFractionOutOfRange(
            generate_pair_command.reference_ancestry_fraction,
        ));
    }

    if generate_pair_command.gap_length_mean.is_nan() {
        return Err(Error::GapLengthMeanIsNaN);
    }
    if generate_pair_command.gap_length_mean < 1.0
        || generate_pair_command.gap_length_mean > generate_pair_command.ancestor_length as f64
    {
        return Err(Error::GapLengthMeanOutOfRange {
            actual: generate_pair_command.gap_length_mean,
            minimum: 1.0,
            maximum: generate_pair_command.ancestor_length as f64,
        });
    }

    let mut input = BufReader::new(File::open(&generate_pair_command.model)?);
    let n: usize = ciborium::from_reader(&mut input)?;
    let alphabet: String = ciborium::from_reader(&mut input)?;
    let alphabet = CliAlphabet::from_str(&alphabet, false).map_err(Error::UnsupportedAlphabet)?;

    if generate_pair_command.ancestor_length < n {
        return Err(Error::LengthLowerThanN {
            length: generate_pair_command.ancestor_length,
            n,
        });
    }

    call::<GeneratePair>(alphabet, n, (input, generate_pair_command))
}

struct GeneratePair;

impl ChooseAlphabetAndN for GeneratePair {
    type Arguments = (BufReader<File>, GeneratePairCommand);

    type Return = ();

    fn call<
        const N: usize,
        const ALPHABET_SIZE: usize,
        BitArrayType: BitViewSized
            + BitStore
            + BitView<Store = BitArrayType>
            + Serialize
            + for<'de> Deserialize<'de>,
        AlphabetType: 'static + Alphabet + IntoCliAlphabet,
    >(
        (input, generate_pair_command): Self::Arguments,
    ) -> Result<Self::Return>
    where
        [u32; ALPHABET_SIZE]: Serialize + for<'de> Deserialize<'de>,
    {
        // Load model.
        let model: NGramModel<N, ALPHABET_SIZE, AlphabetType, BitArrayType> =
            ciborium::from_reader(input)?;

        // Initialise random number generator.
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(generate_pair_command.random_seed);

        // Generate ancestor.
        let ancestor: DefaultGenome<_> =
            model.generate_sequence(generate_pair_command.ancestor_length, &mut rng)?;
        #[expect(unused)]
        let ancestor = if let Some(ancestor_output) = &generate_pair_command.ancestor_output {
            let records = [FastaRecord {
                id: "ancestor".to_string(),
                comment: String::new(),
                sequence_handle: ancestor,
            }];
            write_fasta_file(ancestor_output, &records, &HandleSequenceStore::new())?;
            let [ancestor] = records;
            ancestor.sequence_handle
        } else {
            ancestor
        };

        todo!()
    }
}
