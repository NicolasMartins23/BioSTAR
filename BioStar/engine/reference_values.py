from __future__ import annotations

from BioStar.engine.biochemistry import AminoAcidData


# Immutable biochemical reference data bundled directly with the engine.
# The engine has no runtime database dependency.
AMINO_ACIDS: dict[str, AminoAcidData] = {
    "A": AminoAcidData("A", 71.08, 0.125, 0.710, 0.415, None, None, None),
    "C": AminoAcidData("C", 103.14, -0.815, 0.035, 0.595, 8.18, None, 8.18),
    "D": AminoAcidData("D", 115.09, -1.695, 0.505, 0.485, 3.65, None, 3.65),
    "E": AminoAcidData("E", 129.12, -1.505, 0.755, 0.315, 4.25, None, 4.25),
    "F": AminoAcidData("F", 147.17, 1.925, 0.565, 0.715, None, None, None),
    "G": AminoAcidData("G", 57.05, -0.175, 0.285, 0.625, None, None, None),
    "H": AminoAcidData("H", 137.14, -1.015, 0.500, 0.435, 6.00, None, 6.00),
    "I": AminoAcidData("I", 113.16, -0.225, 0.540, 0.800, None, None, None),
    "K": AminoAcidData("K", 128.18, -2.855, 0.580, 0.370, None, None, 10.53),
    "L": AminoAcidData("L", 113.16, 1.095, 0.605, 0.650, None, None, None),
    "M": AminoAcidData("M", 131.19, 1.765, 0.725, 0.525, None, None, None),
    "N": AminoAcidData("N", 114.10, -1.285, 0.335, 0.445, None, None, None),
    "P": AminoAcidData("P", 97.12, -0.785, 0.285, 0.275, None, None, None),
    "Q": AminoAcidData("Q", 128.13, -1.285, 0.555, 0.550, None, None, None),
    "R": AminoAcidData("R", 156.19, -1.655, 0.490, 0.465, None, None, 12.48),
    "S": AminoAcidData("S", 87.08, -0.095, 0.385, 0.375, None, None, None),
    "T": AminoAcidData("T", 101.11, -0.675, 0.415, 0.595, None, None, None),
    "V": AminoAcidData("V", 99.13, -0.125, 0.530, 0.850, None, None, None),
    "W": AminoAcidData("W", 186.21, -1.865, 0.540, 0.685, None, None, None),
    "Y": AminoAcidData("Y", 163.17, 1.345, 0.345, 0.735, 10.07, None, 10.07),
}

AROMATIC = frozenset({"F", "W", "Y"})
NONPOLAR = frozenset({"A", "G", "I", "L", "M", "P", "V"})
POLAR = frozenset(set(AMINO_ACIDS) - AROMATIC - NONPOLAR)
POSITIVE = frozenset({"H", "K", "R"})
NEGATIVE = frozenset({"D", "E"})

DNA_CODONS: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

RNA_CODONS = {codon.replace("T", "U"): amino_acid for codon, amino_acid in DNA_CODONS.items()}
DNA_STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
RNA_STOP_CODONS = frozenset({"UAA", "UAG", "UGA"})
DNA_START_CODON = "ATG"
RNA_START_CODON = "AUG"

WATER_MASS = 18.015
N_TERM_PKA = 8.00
C_TERM_PKA = 3.10
