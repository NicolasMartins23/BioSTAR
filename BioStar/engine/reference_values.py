from __future__ import annotations

from BioStar.engine.biochemistry import AminoAcidData


# Immutable biochemical reference data bundled directly with the engine.
# Values mirror the legacy BioSTAR biochemical tables; no runtime database is used.
AMINO_ACIDS: dict[str, AminoAcidData] = {
    "A": AminoAcidData("A", 89.09, 1.8, 1.45, 0.97, 2.35, 9.87, None),
    "C": AminoAcidData("C", 121.15, 2.5, 0.77, 1.30, 1.96, 10.28, None),
    "D": AminoAcidData("D", 133.10, -3.5, 1.01, 0.54, 1.88, 9.60, 3.65),
    "E": AminoAcidData("E", 147.13, -3.5, 1.53, 0.37, 2.19, 9.67, 4.25),
    "F": AminoAcidData("F", 165.19, 2.8, 1.13, 1.38, 2.58, 9.24, None),
    "G": AminoAcidData("G", 75.07, -0.4, 0.57, 0.75, 2.34, 9.60, None),
    "H": AminoAcidData("H", 155.16, -3.2, 1.24, 0.87, 1.80, 9.33, 6.04),
    "I": AminoAcidData("I", 131.18, 4.5, 1.00, 1.60, 2.36, 9.60, None),
    "K": AminoAcidData("K", 146.19, -3.9, 1.07, 0.74, 2.18, 9.60, 10.53),
    "L": AminoAcidData("L", 131.18, 3.8, 1.34, 1.22, 2.36, 9.60, None),
    "M": AminoAcidData("M", 149.21, 1.9, 1.20, 1.05, 2.28, 9.21, None),
    "N": AminoAcidData("N", 132.12, -3.5, 0.73, 0.65, 2.02, 8.80, None),
    "P": AminoAcidData("P", 115.13, -1.6, 0.59, 0.62, 1.99, 10.60, None),
    "Q": AminoAcidData("Q", 146.15, -3.5, 1.17, 1.00, 2.17, 9.13, None),
    "R": AminoAcidData("R", 174.20, -4.5, 0.79, 0.90, 2.17, 9.04, 12.48),
    "S": AminoAcidData("S", 105.09, -0.8, 0.82, 0.75, 2.21, 9.15, None),
    "T": AminoAcidData("T", 119.12, -0.7, 0.83, 1.19, 2.09, 9.10, None),
    "V": AminoAcidData("V", 117.15, 4.2, 1.06, 1.70, 2.32, 9.62, None),
    "W": AminoAcidData("W", 204.23, -0.9, 1.08, 1.37, 2.38, 9.39, None),
    "Y": AminoAcidData("Y", 181.19, -1.3, 0.69, 1.47, 2.20, 9.11, None),
}

AROMATIC = frozenset({"F", "W", "Y"})
NONPOLAR = frozenset({"A", "C", "G", "I", "L", "M", "P", "V"})
POLAR = frozenset({"D", "E", "H", "K", "N", "Q", "R", "S", "T"})
POSITIVE = frozenset({"K", "R", "H"})
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

WATER_MASS = 18.01528
N_TERM_PKA = 7.7
C_TERM_PKA = 3.5
