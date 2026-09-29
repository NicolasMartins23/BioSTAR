from BioStar import DNA, RNA


def test_dna_uses_bundled_codon_table() -> None:
    dna = DNA("ATGGGC")

    assert dna.get_peptide_sequence() == "MG"
    assert dna.rna_sequence() == "AUGGGC"


def test_rna_uses_bundled_codon_table() -> None:
    rna = RNA("AUGGGC")

    assert rna.get_peptide_sequence() == "MG"
    assert rna.dna_sequence() == "ATGGGC"


def test_dna_and_rna_use_bundled_stop_codons() -> None:
    dna = DNA("ATGTAA")
    rna = RNA("AUGUAA")

    assert dna.get_peptide_sequence() == "M"
    assert rna.trim_on_stop_codon() == "AUGUAA"
