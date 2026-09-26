from BioStar.analysis import CompareNucleotideSequence


def test_sequence_comparison() -> None:
    comparison: CompareNucleotideSequence = CompareNucleotideSequence("ATG", "ATA")
    result: list[dict[str, object]] = comparison.compare()
    assert len(result) == 2
    assert result[0]["mutation_type"] == "missense"
    assert result[1]["mutation_type"] == "no_mutation"
