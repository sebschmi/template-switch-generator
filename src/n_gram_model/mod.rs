use std::collections::BTreeMap;

use compact_genome::{
    implementation::bit_array_kmer::{BitArrayKmer, BitStore, BitViewSized},
    interface::{
        alphabet::{Alphabet, AlphabetCharacter},
        sequence::GenomeSequence,
    },
};

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
        'item,
        SubsequenceType: 'item + GenomeSequence<AlphabetType, SubsequenceType> + ?Sized,
    >(
        sequences: impl IntoIterator<Item = &'item SubsequenceType>,
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
                    abundances[successor.index()] += 1;
                } else {
                    let mut abundances = [0; ALPHABET_SIZE];
                    abundances[successor.index()] = 1;
                    result.model.insert(kmer, abundances);
                }
            }
        }

        result
    }
}
