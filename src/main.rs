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
use log::{info, LevelFilter};
use n_gram_model::NGramModel;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use sequence_modifier::{
    template_switch_overlap_detector::TemplateSwitchOverlapDetector, SequenceModifier,
    SequenceModifierPair,
};
use serde::{Deserialize, Serialize};
use simplelog::{ColorChoice, TermLogger, TerminalMode};

mod choose_alphabet_and_n;
mod cli;
mod error;
mod n_gram_model;
mod sequence_modifier;

fn main() {
    TermLogger::init(
        LevelFilter::Info,
        Default::default(),
        TerminalMode::Mixed,
        ColorChoice::Auto,
    )
    .unwrap();

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
        info!("Loading sequences...");
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
        info!("Creating model...");
        let model = NGramModel::<N, ALPHABET_SIZE, _, BitArrayType>::from_sequences(sequences);

        // Write model parameters and model.
        info!("Storing model...");
        let mut output = BufWriter::new(File::create(&create_model_command.output)?);
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
    generate_pair_command.verify()?;

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

        // Derive reference and query from ancestor.
        let mut reference = ancestor.clone();
        let mut query = ancestor.clone();

        let SequenceModifierPair {
            mut reference_modifier,
            mut query_modifier,
        } = SequenceModifier::new_modifier_pair(
            generate_pair_command.reference_ancestry_fraction,
            generate_pair_command.sequence_modification_amount,
            generate_pair_command.sequence_modification_parameters,
            &mut rng,
        );

        let mut template_switch_overlap_detector = TemplateSwitchOverlapDetector::new(
            &generate_pair_command.sequence_modification_parameters,
        );
        reference_modifier.apply(
            &mut reference,
            &mut template_switch_overlap_detector,
            &mut rng,
        )?;
        template_switch_overlap_detector.clear_modification_stack();
        query_modifier.apply(&mut query, &mut template_switch_overlap_detector, &mut rng)?;

        // Write sequences.
        write_fasta_file(
            &generate_pair_command.output,
            &[
                FastaRecord {
                    id: "reference".to_string(),
                    comment: String::new(),
                    sequence_handle: reference,
                },
                FastaRecord {
                    id: "query".to_string(),
                    comment: String::new(),
                    sequence_handle: query,
                },
            ],
            &HandleSequenceStore::new(),
        )?;

        Ok(())
    }
}
