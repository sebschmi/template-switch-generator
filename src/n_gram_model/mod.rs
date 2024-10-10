use std::collections::BTreeMap;

use compact_genome::{
    implementation::bit_array_kmer::{BitArrayKmer, BitStore, BitView, BitViewSized},
    interface::{
        alphabet::{Alphabet, AlphabetCharacter},
        k_mer::OwnedKmer,
        sequence::{GenomeSequence, OwnedGenomeSequence},
    },
};
use rand::{
    distributions::{Uniform, WeightedIndex},
    prelude::Distribution,
    Rng,
};

use crate::error::{Error, Result};

mod serde;

pub struct NGramModel<
    const N: usize,
    const ALPHABET_SIZE: usize,
    AlphabetType: Alphabet,
    BitArrayType: BitViewSized + BitStore,
> {
    model: BTreeMap<BitArrayKmer<N, AlphabetType, BitArrayType>, [u32; ALPHABET_SIZE]>,
}

impl<
        const N: usize,
        const ALPHABET_SIZE: usize,
        AlphabetType: Alphabet,
        BitArrayType: BitViewSized + BitStore,
    > NGramModel<N, ALPHABET_SIZE, AlphabetType, BitArrayType>
{
    pub fn from_sequences<
        SequenceType: GenomeSequence<AlphabetType, SubsequenceType>,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        sequences: impl IntoIterator<Item = SequenceType>,
    ) -> Self {
        assert_eq!(ALPHABET_SIZE, AlphabetType::SIZE);
        let mut result = Self {
            model: Default::default(),
        };

        for sequence in sequences {
            for offset in 0..sequence.len() - N - 1 {
                let kmer = BitArrayKmer::from_iter(sequence[offset..offset + N].iter().cloned());
                let successor = sequence[offset + N].clone();

                if let Some(abundances) = result.model.get_mut(&kmer) {
                    abundances[successor.index()] =
                        abundances[successor.index()].checked_add(1).unwrap();
                } else {
                    let mut abundances = [0; ALPHABET_SIZE];
                    abundances[successor.index()] = 1;
                    result.model.insert(kmer, abundances);
                }
            }
        }

        result
    }

    pub fn generate_sequence<
        SequenceType: OwnedGenomeSequence<AlphabetType, SubsequenceType>,
        SubsequenceType: GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        &self,
        length: usize,
        rng: &mut impl Rng,
    ) -> Result<SequenceType>
    where
        BitArrayType: BitView<Store = BitArrayType>,
    {
        if length < N {
            return Err(Error::LengthLowerThanN { length, n: N });
        }

        let kmer_sampler = WeightedIndex::<usize>::new(self.model.values().map(|abundances| {
            abundances
                .iter()
                .map(|abundance| *abundance as usize)
                .sum::<usize>()
        }))
        .map_err(|_| Error::EmptyModel)?;
        let generator = NGramSequenceGenerator::new(self, rng, kmer_sampler);
        Ok(SequenceType::from_iter(generator.take(length)))
    }
}

struct NGramSequenceGenerator<
    'model,
    'rng,
    const N: usize,
    const ALPHABET_SIZE: usize,
    AlphabetType: Alphabet,
    BitArrayType: BitViewSized + BitStore + BitView<Store = BitArrayType>,
    RandomNumberGenerator: Rng,
> {
    kmer: Option<BitArrayKmer<N, AlphabetType, BitArrayType>>,
    next_index: usize,
    model: &'model NGramModel<N, ALPHABET_SIZE, AlphabetType, BitArrayType>,
    rng: &'rng mut RandomNumberGenerator,
    kmer_sampler: WeightedIndex<usize>,
}

impl<
        'model,
        'rng,
        const N: usize,
        const ALPHABET_SIZE: usize,
        AlphabetType: Alphabet,
        BitArrayType: BitViewSized + BitStore + BitView<Store = BitArrayType>,
        RandomNumberGenerator: Rng,
    >
    NGramSequenceGenerator<
        'model,
        'rng,
        N,
        ALPHABET_SIZE,
        AlphabetType,
        BitArrayType,
        RandomNumberGenerator,
    >
{
    fn new(
        model: &'model NGramModel<N, ALPHABET_SIZE, AlphabetType, BitArrayType>,
        rng: &'rng mut RandomNumberGenerator,
        kmer_sampler: WeightedIndex<usize>,
    ) -> Self {
        Self {
            kmer: None,
            next_index: 0,
            model,
            rng,
            kmer_sampler,
        }
    }
}

impl<
        const N: usize,
        const ALPHABET_SIZE: usize,
        AlphabetType: Alphabet,
        BitArrayType: BitViewSized + BitStore + BitView<Store = BitArrayType>,
        RandomNumberGenerator: Rng,
    > Iterator
    for NGramSequenceGenerator<
        '_,
        '_,
        N,
        ALPHABET_SIZE,
        AlphabetType,
        BitArrayType,
        RandomNumberGenerator,
    >
{
    type Item = AlphabetType::CharacterType;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(kmer) = &mut self.kmer {
            if self.next_index < N {
                let result = kmer[self.next_index].clone();
                self.next_index += 1;
                Some(result)
            } else if let Some(abundances) = self.model.model.get(kmer) {
                let sum: u32 = abundances.iter().cloned().sum();
                let distribution = Uniform::new(0, sum);
                let sample = distribution.sample(self.rng);

                let mut index = usize::MAX;
                let mut current_sum = 0;
                for (current_index, value) in abundances.iter().cloned().enumerate() {
                    current_sum += value;
                    if sample < current_sum {
                        index = current_index;
                        break;
                    }
                }
                debug_assert_ne!(index, usize::MAX);

                let character = AlphabetType::CharacterType::from_index(index).unwrap();
                self.kmer = Some(kmer.successor(character.clone()));
                Some(character)
            } else {
                self.kmer = None;
                self.next()
            }
        } else {
            self.kmer = Some(
                self.model
                    .model
                    .keys()
                    .nth(self.kmer_sampler.sample(self.rng))
                    .unwrap()
                    .clone(),
            );
            self.next_index = 0;
            self.next()
        }
    }
}
