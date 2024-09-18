use compact_genome::{
    implementation::bit_array_kmer::{BitStore, BitViewSized},
    interface::alphabet::Alphabet,
};
use serde::{Deserialize, Serialize};

use super::NGramModel;

impl<
        const N: usize,
        const ALPHABET_SIZE: usize,
        AlphabetType: Alphabet,
        BitArrayType: BitViewSized + BitStore + Serialize,
    > Serialize for NGramModel<N, ALPHABET_SIZE, AlphabetType, BitArrayType>
where
    [u32; ALPHABET_SIZE]: Serialize,
{
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        self.model.serialize(serializer)
    }
}

impl<
        'de,
        const N: usize,
        const ALPHABET_SIZE: usize,
        AlphabetType: Alphabet,
        BitArrayType: BitViewSized + BitStore + Deserialize<'de>,
    > Deserialize<'de> for NGramModel<N, ALPHABET_SIZE, AlphabetType, BitArrayType>
where
    [u32; ALPHABET_SIZE]: Deserialize<'de>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        Ok(Self {
            model: Deserialize::deserialize(deserializer)?,
        })
    }
}
