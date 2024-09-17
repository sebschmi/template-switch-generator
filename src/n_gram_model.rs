use std::collections::BTreeMap;

use compact_genome::interface::alphabet::Alphabet;

#[expect(unused)]
pub struct NGramModel<const N: usize, AlphabetType: Alphabet> {
    model: BTreeMap<NGramEntry<N, AlphabetType>, u32>,
}

#[expect(unused)]
struct NGramEntry<const N: usize, AlphabetType: Alphabet> {
    identifier: [AlphabetType::CharacterType; N],
    successor: AlphabetType::CharacterType,
}
