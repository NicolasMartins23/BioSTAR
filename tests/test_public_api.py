from BioStar import DNA, Protein, RNA


def test_public_imports() -> None:
    assert DNA("ATG").sequence == "ATG"
    assert RNA("AUG").sequence == "AUG"
    assert Protein("ACDE").sequence == "ACDE"
