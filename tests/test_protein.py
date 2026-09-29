from BioStar import Protein


def test_protein_uses_bundled_reference_data() -> None:
    protein = Protein("MG")

    assert protein.sequence == "MG"
    assert protein.molecular_weight() == 206.26


def test_protein_analysis_methods_use_reference_data() -> None:
    protein = Protein("ACDEFGHIKLMNPQRSTVWY")

    assert protein.sequence_size == 20
    assert protein.aromacity() == 0.15
    assert protein.hydrophobic_index() == -0.49
