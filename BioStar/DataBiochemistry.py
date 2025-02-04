CODON_SIZE: int = 3
WATER_MASS: float = 18.01528
N_TERM_PKA: float = 7.7
C_TERM_PKA: float = 3.5

AMINOACIDS: list[str] = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
AMINOACIDS_AROMATIC: list[str] = ["F", "W", "Y"]
AMINOACIDS_NONPOLAR: list[str] = ["A", "C", "G", "I", "L", "M", "P", "V"]
AMINOACIDS_POLAR: list[str] = ["D", "E", "H", "K", "N", "Q", "R", "S", "T", "Q"]
AMINOACIDS_POSITIVE: list[str] = ["K", "R", "H"]  # Positively charged amino acids
AMINOACIDS_NEGATIVE: list[str] = ["D", "E"]  # Negatively charged amino acids

AMINOACID_TABLE_TO_CODON_DNA: dict = {
    '*': ['TAA', 'TGA', 'TAG'],
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAC', 'GAT'],
    'E': ['GAA', 'GAG'],
    'F': ['TTC', 'TTT'],
    'G': ['GGA', 'GGC', 'GGG', 'GGT'],
    'H': ['CAC', 'CAT'],
    'I': ['ATA', 'ATC', 'ATT'],
    'K': ['AAA', 'AAG'],
    'L': ['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'],
    'M': ['ATG'],
    'N': ['AAC', 'AAT'],
    'P': ['CCA', 'CCC', 'CCG', 'CCT'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
    'S': ['AGC', 'AGT', 'TCA', 'TCC', 'TCG', 'TCT'],
    'T': ['ACA', 'ACC', 'ACG', 'ACT'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'W': ['TGG'],
    'Y': ['TAC', 'TAT']
}

AMINOACID_TABLE_TO_CODON_RNA: dict = {
    '*': ['UAA', 'UGA', 'UAG'],
    'A': ['GCA', 'GCC', 'GCG', 'GCU'],
    'C': ['UGC', 'UGU'],
    'D': ['GAC', 'GAU'],
    'E': ['GAA', 'GAG'],
    'F': ['UUC', 'UUU'],
    'G': ['GGA', 'GGC', 'GGG', 'GGU'],
    'H': ['CAC', 'CAU'],
    'I': ['AUA', 'AUC', 'AUU'],
    'K': ['AAA', 'AAG'],
    'L': ['CUA', 'CUC', 'CUG', 'CUU', 'UUA', 'UUG'],
    'M': ['AUG'],
    'N': ['AAC', 'AAU'],
    'P': ['CCA', 'CCC', 'CCG', 'CCU'],
    'Q': ['CAA', 'CAG'],
    'R': ['AGA', 'AGG', 'CGA', 'CGC', 'CGG', 'CGU'],
    'S': ['AGC', 'AGU', 'UCA', 'UCC', 'UCG', 'UCU'],
    'T': ['ACA', 'ACC', 'ACG', 'ACU'],
    'V': ['GUA', 'GUC', 'GUG', 'GUU'],
    'W': ['UGG'],
    'Y': ['UAC', 'UAU']
}

TABLE_DNA_CODON_TO_AMINOACID: dict = {
#A
# A           C           G           T
'AAA': 'K', 'ACA': 'T', 'AGA': 'R', 'ATA': 'I', #A
'AAC': 'N', 'ACC': 'T', 'AGC': 'S', 'ATC': 'I', #C
'AAG': 'K', 'ACG': 'T', 'AGG': 'R', 'ATG': 'M', #G
'AAT': 'N', 'ACT': 'T', 'AGT': 'S', 'ATT': 'I', #T

#C
# A           C           G           T
'CAA': 'Q', 'CCA': 'P', 'CGA': 'R', 'CTA': 'L', #A
'CAC': 'H', 'CCC': 'P', 'CGC': 'R', 'CTC': 'L', #C
'CAG': 'Q', 'CCG': 'P', 'CGG': 'R', 'CTG': 'L', #G
'CAT': 'H', 'CCT': 'P', 'CGT': 'R', 'CTT': 'L', #T

#G
# A           C           G           T
'GAA': 'E', 'GCA': 'A', 'GGA': 'G', 'GTA': 'V', #A
'GAC': 'D', 'GCC': 'A', 'GGC': 'G', 'GTC': 'V', #C
'GAG': 'E', 'GCG': 'A', 'GGG': 'G', 'GTG': 'V', #G
'GAT': 'D', 'GCT': 'A', 'GGT': 'G', 'GTT': 'V', #T

#T
# A           C           G           T
'TAA': '*', 'TCA': 'S', 'TGA': '*', 'TTA': 'L', #A
'TAC': 'Y', 'TCC': 'S', 'TGC': 'C', 'TTC': 'F', #C
'TAG': '*', 'TCG': 'S', 'TGG': 'W', 'TTG': 'L', #G
'TAT': 'Y', 'TCT': 'S', 'TGT': 'C', 'TTT': 'F', #T
}


TABLE_RNA_CODON_TO_AMINOACID: dict = {
#A
# A           C           G           U
'AAA': 'K', 'ACA': 'T', 'AGA': 'R', 'AUA': 'I', #A
'AAC': 'N', 'ACC': 'T', 'AGC': 'S', 'AUC': 'I', #C
'AAG': 'K', 'ACG': 'T', 'AGG': 'R', 'AUG': 'M', #G
'AAU': 'N', 'ACU': 'T', 'AGU': 'S', 'AUU': 'I', #U

#C
# A           C           G           U
'CAA': 'Q', 'CCA': 'P', 'CGA': 'R', 'CUA': 'L', #A
'CAC': 'H', 'CCC': 'P', 'CGC': 'R', 'CUC': 'L', #C
'CAG': 'Q', 'CCG': 'P', 'CGG': 'R', 'CUG': 'L', #G
'CAU': 'H', 'CCU': 'P', 'CGU': 'R', 'CUU': 'L', #U

#G
# A           C           G           U
'GAA': 'E', 'GCA': 'A', 'GGA': 'G', 'GUA': 'V', #A
'GAC': 'D', 'GCC': 'A', 'GGC': 'G', 'GUC': 'V', #C
'GAG': 'E', 'GCG': 'A', 'GGG': 'G', 'GUG': 'V', #G
'GAU': 'D', 'GCU': 'A', 'GGU': 'G', 'GUU': 'V', #U

#U
# A           C           G           U
'UAA': '_', 'UCA': 'S', 'UGA': '_', 'UUA': 'L', #A
'UAC': 'Y', 'UCC': 'S', 'UGC': 'C', 'UUC': 'F', #C
'UAG': '_', 'UCG': 'S', 'UGG': 'W', 'UUG': 'L', #G
'UAU': 'Y', 'UCU': 'S', 'UGU': 'C', 'UUU': 'F', #U
}

START_CODON_DNA = "ATG"
STOP_CODON_DNA = ["TAA", "TAG", "TGA"]

START_CODON_TNA = "ATG"
STOP_CODON_RNA = ["UAA", "UAG", "UGA"]

AMINOACID_TABLE = {
    "A": {
        "single_letter": "A",
        "abbreviation": "Ala",
        "name": "Alanine",
        "dna_codons": ["GCA", "GCC", "GCG", "GCT"],
        "rna_codons": ["GCA", "GCC", "GCG", "GCU"],
        "weight": 89.09,
        "hydrophobicity": 1.8,
        "alpha_helix": 1.45,
        "beta_sheet": 0.97
    },
    "C": {
        "single_letter": "C",
        "abbreviation": "Cys",
        "name": "Cysteine",
        "dna_codons": ["TGC", "TGT"],
        "rna_codons": ["UGC", "UGU"],
        "weight": 121.15,
        "hydrophobicity": 2.5,
        "alpha_helix": 0.77,
        "beta_sheet": 1.30
    },
    "D": {
        "single_letter": "D",
        "abbreviation": "Asp",
        "name": "Aspartic acid",
        "dna_codons": ["GAC", "GAT"],
        "rna_codons": ["GAC", "GAU"],
        "weight": 133.10,
        "hydrophobicity": -3.5,
        "alpha_helix": 1.01,
        "beta_sheet": 0.54
    },
    "E": {
        "single_letter": "E",
        "abbreviation": "Glu",
        "name": "Glutamic acid",
        "dna_codons": ["GAA", "GAG"],
        "rna_codons": ["GAA", "GAG"],
        "weight": 147.13,
        "hydrophobicity": -3.5,
        "alpha_helix": 1.53,
        "beta_sheet": 0.37
    },
    "F": {
        "single_letter": "F",
        "abbreviation": "Phe",
        "name": "Phenylalanine",
        "dna_codons": ["TTC", "TTT"],
        "rna_codons": ["UUC", "UUU"],
        "weight": 165.19,
        "hydrophobicity": 2.8,
        "alpha_helix": 1.13,
        "beta_sheet": 1.38
    },
    "G": {
        "single_letter": "G",
        "abbreviation": "Gly",
        "name": "Glycine",
        "dna_codons": ["GGA", "GGC", "GGG", "GGT"],
        "rna_codons": ["GGA", "GGC", "GGG", "GGU"],
        "weight": 75.07,
        "hydrophobicity": -0.4,
        "alpha_helix": 0.57,
        "beta_sheet": 0.75
    },
    "H": {
        "single_letter": "H",
        "abbreviation": "His",
        "name": "Histidine",
        "dna_codons": ["CAC", "CAT"],
        "rna_codons": ["CAC", "CAU"],
        "weight": 155.16,
        "hydrophobicity": -3.2,
        "alpha_helix": 1.24,
        "beta_sheet": 0.87
    },
    "I": {
        "single_letter": "I",
        "abbreviation": "Ile",
        "name": "Isoleucine",
        "dna_codons": ["ATA", "ATC", "ATT"],
        "rna_codons": ["AUA", "AUC", "AUU"],
        "weight": 131.18,
        "hydrophobicity": 4.5,
        "alpha_helix": 1.00,
        "beta_sheet": 1.60
    },
    "K": {
        "single_letter": "K",
        "abbreviation": "Lys",
        "name": "Lysine",
        "dna_codons": ["AAA", "AAG"],
        "rna_codons": ["AAA", "AAG"],
        "weight": 146.19,
        "hydrophobicity": -3.9,
        "alpha_helix": 1.07,
        "beta_sheet": 0.74
    },
    "L": {
        "single_letter": "L",
        "abbreviation": "Leu",
        "name": "Leucine",
        "dna_codons": ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
        "rna_codons": ["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"],
        "weight": 131.18,
        "hydrophobicity": 3.8,
        "alpha_helix": 1.34,
        "beta_sheet": 1.22
    },
    "M": {
        "single_letter": "M",
        "abbreviation": "Met",
        "name": "Methionine",
        "dna_codons": ["ATG"],
        "rna_codons": ["AUG"],
        "weight": 149.21,
        "hydrophobicity": 1.9,
        "alpha_helix": 1.20,
        "beta_sheet": 1.05
    },
    "N": {
        "single_letter": "N",
        "abbreviation": "Asn",
        "name": "Asparagine",
        "dna_codons": ["AAC", "AAT"],
        "rna_codons": ["AAC", "AAU"],
        "weight": 132.12,
        "hydrophobicity": -3.5,
        "alpha_helix": 0.73,
        "beta_sheet": 0.65
    },
    "P": {
        "single_letter": "P",
        "abbreviation": "Pro",
        "name": "Proline",
        "dna_codons": ["CCA", "CCC", "CCG", "CCT"],
        "rna_codons": ["CCA", "CCC", "CCG", "CCU"],
        "weight": 115.13,
        "hydrophobicity": -1.6,
        "alpha_helix": 0.59,
        "beta_sheet": 0.62
    },
    "Q": {
        "single_letter": "Q",
        "abbreviation": "Gln",
        "name": "Glutamine",
        "dna_codons": ["CAA", "CAG"],
        "rna_codons": ["CAA", "CAG"],
        "weight": 146.15,
        "hydrophobicity": -3.5,
        "alpha_helix": 1.17,
        "beta_sheet": 1.00
    },
    "R": {
        "single_letter": "R",
        "abbreviation": "Arg",
        "name": "Arginine",
        "dna_codons": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
        "rna_codons": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"],
        "weight": 174.20,
        "hydrophobicity": -4.5,
        "alpha_helix": 0.79,
        "beta_sheet": 0.90
    },
    "S": {
        "single_letter": "S",
        "abbreviation": "Ser",
        "name": "Serine",
        "dna_codons": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
        "rna_codons": ["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"],
        "weight": 105.09,
        "hydrophobicity": -0.8,
        "alpha_helix": 0.82,
        "beta_sheet": 0.75
    },
    "T": {
        "single_letter": "T",
        "abbreviation": "Thr",
        "name": "Threonine",
        "dna_codons": ["ACA", "ACC", "ACG", "ACT"],
        "rna_codons": ["ACA", "ACC", "ACG", "ACU"],
        "weight": 119.12,
        "hydrophobicity": -0.7,
        "alpha_helix": 0.83,
        "beta_sheet": 1.19
    },
    "V": {
        "single_letter": "V",
        "abbreviation": "Val",
        "name": "Valine",
        "dna_codons": ["GTA", "GTC", "GTG", "GTT"],
        "rna_codons": ["GUA", "GUC", "GUG", "GUU"],
        "weight": 117.15,
        "hydrophobicity": 4.2,
        "alpha_helix": 1.06,
        "beta_sheet": 1.70
    },
    "W": {
        "single_letter": "W",
        "abbreviation": "Trp",
        "name": "Tryptophan",
        "dna_codons": ["TGG"],
        "rna_codons": ["UGG"],
        "weight": 204.23,
        "hydrophobicity": -0.9,
        "alpha_helix": 1.08,
        "beta_sheet": 1.37
    },
    "Y": {
        "single_letter": "Y",
        "abbreviation": "Tyr",
        "name": "Tyrosine",
        "dna_codons": ["TAC", "TAT"],
        "rna_codons": ["UAC", "UAU"],
        "weight": 181.19,
        "hydrophobicity": -1.3,
        "alpha_helix": 0.69,
        "beta_sheet": 1.47
    },
    "*": {
        "single_letter": "*",
        "abbreviation": "Stop",
        "name": "Stop Codon",
        "dna_codons": ["TAA", "TAG", "TGA"],
        "rna_codons": ["UAA", "UAG", "UGA"],
        "weight": 0.0,
        "hydrophobicity": 0.0,
        "alpha_helix": 0.0,
        "beta_sheet": 0.0
    }
}
