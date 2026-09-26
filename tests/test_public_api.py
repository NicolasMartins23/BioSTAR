from BioStar import DNA, Protein, RNA


def test_public_imports() -> None:
    assert DNA("atg").sequence == "ATG"
    assert RNA("aug").sequence == "AUG"
    assert Protein("acde").sequence == "ACDE"


def test_dna_transforms() -> None:
    dna: DNA = DNA("ATGTAA")
    assert dna.to_rna().sequence == "AUGUAA"
    assert dna.get_peptide_sequence() == "M"
    assert dna.to_protein().sequence == "M"


def test_protein_empty_sequence_is_safe() -> None:
    protein: Protein = Protein()
    assert protein.molecular_weight() == 0.0
    assert protein.hydrophobic_index() == 0.0
    assert protein.secondary_structure_propensity()["coil"] == 0.0
