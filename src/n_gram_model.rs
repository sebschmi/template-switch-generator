use compact_genome::interface::alphabet::Alphabet;

#[expect(unused)]
pub struct NGramModel<const N: usize, AlphabetType: Alphabet> {
    model: Vec<NGramEntry<N, AlphabetType>>,
}

#[expect(unused)]
struct NGramEntry<const N: usize, AlphabetType: Alphabet> {
    identifier: [AlphabetType::CharacterType; N],
    successor: AlphabetType::CharacterType,
    amount: u32,
}
