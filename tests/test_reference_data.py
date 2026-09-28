from BioStar.engine import get_default_biochemistry


def test_biochemical_reference_data_loads_from_sqlite() -> None:
    data = get_default_biochemistry()

    assert len(data.amino_acids) == 20
    assert len(data.dna_codons) == 64
    assert len(data.rna_codons) == 64
    assert data.dna_start_codon == "ATG"
    assert data.rna_start_codon == "AUG"
    assert data.dna_stop_codons == frozenset({"TAA", "TAG", "TGA"})
    assert data.rna_stop_codons == frozenset({"UAA", "UAG", "UGA"})
