from BioStar.engine import get_reference_data


def test_biochemical_reference_data_loads_from_sqlite() -> None:
    data = get_reference_data()

    assert len(data.amino_acid_symbols()) == 20
    assert len(data.codon_table("DNA")) == 64
    assert len(data.codon_table("RNA")) == 64
    assert data.start_codon("DNA") == "ATG"
    assert data.start_codon("RNA") == "AUG"
    assert data.stop_codons("DNA") == frozenset({"TAA", "TAG", "TGA"})
    assert data.stop_codons("RNA") == frozenset({"UAA", "UAG", "UGA"})
