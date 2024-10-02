use compact_genome::{
    implementation::{
        alphabets::dna_alphabet::DnaAlphabet,
        bit_array_kmer::{BitStore, BitView, BitViewSized},
    },
    interface::alphabet::Alphabet,
};
use serde::{Deserialize, Serialize};

use crate::{
    cli::{CliAlphabet, IntoCliAlphabet},
    error::{Error, Result},
};

pub fn call<Function: ChooseAlphabetAndN>(
    alphabet: CliAlphabet,
    n: usize,
    arguments: Function::Arguments,
) -> Result<Function::Return> {
    match alphabet {
        CliAlphabet::Dna => {
            with_alphabet::<{ DnaAlphabet::SIZE }, DnaAlphabet, Function>(n, arguments)
        }
    }
}

fn with_alphabet<
    const ALPHABET_SIZE: usize,
    AlphabetType: 'static + Alphabet + IntoCliAlphabet,
    Function: ChooseAlphabetAndN,
>(
    n: usize,
    arguments: Function::Arguments,
) -> Result<Function::Return>
where
    [u32; ALPHABET_SIZE]: Serialize + for<'de> Deserialize<'de>,
{
    match n {
        0 => with_alphabet_and_n::<0, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        1 => with_alphabet_and_n::<1, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        2 => with_alphabet_and_n::<2, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        3 => with_alphabet_and_n::<3, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        4 => with_alphabet_and_n::<4, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        5 => with_alphabet_and_n::<5, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        6 => with_alphabet_and_n::<6, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        7 => with_alphabet_and_n::<7, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        8 => with_alphabet_and_n::<8, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        9 => with_alphabet_and_n::<9, ALPHABET_SIZE, AlphabetType, Function>(arguments),
        n => Err(Error::UnsupportedN(n)),
    }
}

fn with_alphabet_and_n<
    const N: usize,
    const ALPHABET_SIZE: usize,
    AlphabetType: 'static + Alphabet + IntoCliAlphabet,
    Function: ChooseAlphabetAndN,
>(
    arguments: Function::Arguments,
) -> Result<Function::Return>
where
    [u32; ALPHABET_SIZE]: Serialize + for<'de> Deserialize<'de>,
{
    let n_gram_bit_width = (ALPHABET_SIZE + 1).ilog2() as usize;
    let bit_width = n_gram_bit_width * N;

    if bit_width <= 8 {
        Function::call::<N, ALPHABET_SIZE, u8, AlphabetType>(arguments)
    } else if bit_width <= 16 {
        Function::call::<N, ALPHABET_SIZE, u16, AlphabetType>(arguments)
    } else if bit_width <= 32 {
        Function::call::<N, ALPHABET_SIZE, u32, AlphabetType>(arguments)
    } else if bit_width <= 64 {
        Function::call::<N, ALPHABET_SIZE, u64, AlphabetType>(arguments)
    } else {
        Err(Error::UnsupportedN(N))
    }
}

pub trait ChooseAlphabetAndN {
    type Arguments;
    type Return;

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
        arguments: Self::Arguments,
    ) -> Result<Self::Return>
    where
        [u32; ALPHABET_SIZE]: Serialize + for<'de> Deserialize<'de>;
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_ilog2() {
        assert_eq!(1usize.ilog2(), 0);

        assert_eq!(2usize.ilog2(), 1);
        assert_eq!(3usize.ilog2(), 1);

        assert_eq!(4usize.ilog2(), 2);
        assert_eq!(5usize.ilog2(), 2);
        assert_eq!(6usize.ilog2(), 2);
        assert_eq!(7usize.ilog2(), 2);

        assert_eq!(8usize.ilog2(), 3);
        assert_eq!(9usize.ilog2(), 3);
    }
}
