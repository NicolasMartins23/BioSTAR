
AMINOACIDS = "ACDEFGHIKLMNPQRSTVWY"
AMINOACIDS_AROMATIC = "FWY"
AMINOACIDS_NONPOLAR = "ACGILMPV"
AMINOACIDS_POLAR = "DEHKNQRSTQ"

CODON_SIZE = int(3)

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

AMINOACID_WEIGHT = {
    "A": 89.1,
    "R": 174.2,
    "N": 132.1,
    "D": 133.1,
    "C": 121.2,
    "E": 147.1,
    "Q": 146.2,
    "G": 75.1,
    "H": 155.2,
    "I": 131.2,
    "L": 131.2,
    "K": 146.2,
    "M": 149.2,
    "F": 165.2,
    "P": 115.1,
    "S": 105.1,
    "T": 119.1,
    "W": 204.2,
    "Y": 181.2,
    "V": 117.1,
    "*": 0
}


AMINOACID_TABLE = {
    "A": {
        "single_letter": "A",
        "abbreviation": "Ala",
        "name": "Alanine",
        "dna_codons": ["GCA", "GCC", "GCG", "GCT"],
        "rna_codons": ["GCA", "GCC", "GCG", "GCU"]
    },
    "C": {
        "single_letter": "C",
        "abbreviation": "Cys",
        "name": "Cysteine",
        "dna_codons": ["TGC", "TGT"],
        "rna_codons": ["UGC", "UGU"]
    },
    "D": {
        "single_letter": "D",
        "abbreviation": "Asp",
        "name": "Aspartic acid",
        "dna_codons": ["GAC", "GAT"],
        "rna_codons": ["GAC", "GAU"]
    },
    "E": {
        "single_letter": "E",
        "abbreviation": "Glu",
        "name": "Glutamic acid",
        "dna_codons": ["GAA", "GAG"],
        "rna_codons": ["GAA", "GAG"]
    },
    "F": {
        "single_letter": "F",
        "abbreviation": "Phe",
        "name": "Phenylalanine",
        "dna_codons": ["TTC", "TTT"],
        "rna_codons": ["UUC", "UUU"]
    },
    "G": {
        "single_letter": "G",
        "abbreviation": "Gly",
        "name": "Glycine",
        "dna_codons": ["GGA", "GGC", "GGG", "GGT"],
        "rna_codons": ["GGA", "GGC", "GGG", "GGU"]
    },
    "H": {
        "single_letter": "H",
        "abbreviation": "His",
        "name": "Histidine",
        "dna_codons": ["CAC", "CAT"],
        "rna_codons": ["CAC", "CAU"]
    },
    "I": {
        "single_letter": "I",
        "abbreviation": "Ile",
        "name": "Isoleucine",
        "dna_codons": ["ATA", "ATC", "ATT"],
        "rna_codons": ["AUA", "AUC", "AUU"]
    },
    "K": {
        "single_letter": "K",
        "abbreviation": "Lys",
        "name": "Lysine",
        "dna_codons": ["AAA", "AAG"],
        "rna_codons": ["AAA", "AAG"]
    },
    "L": {
        "single_letter": "L",
        "abbreviation": "Leu",
        "name": "Leucine",
        "dna_codons": ["CTA", "CTC", "CTG", "CTT", "TTA", "TTG"],
        "rna_codons": ["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"]
    },
    "M": {
        "single_letter": "M",
        "abbreviation": "Met",
        "name": "Methionine",
        "dna_codons": ["ATG"],
        "rna_codons": ["AUG"]
    },
    "N": {
        "single_letter": "N",
        "abbreviation": "Asn",
        "name": "Asparagine",
        "dna_codons": ["AAC", "AAT"],
        "rna_codons": ["AAC", "AAU"]
    },
    "P": {
        "single_letter": "P",
        "abbreviation": "Pro",
        "name": "Proline",
        "dna_codons": ["CCA", "CCC", "CCG", "CCT"],
        "rna_codons": ["CCA", "CCC", "CCG", "CCU"]
    },
    "Q": {
        "single_letter": "Q",
        "abbreviation": "Gln",
        "name": "Glutamine",
        "dna_codons": ["CAA", "CAG"],
        "rna_codons": ["CAA", "CAG"]
    },
    "R": {
        "single_letter": "R",
        "abbreviation": "Arg",
        "name": "Arginine",
        "dna_codons": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGT"],
        "rna_codons": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"]
    },
    "S": {
        "single_letter": "S",
        "abbreviation": "Ser",
        "name": "Serine",
        "dna_codons": ["AGC", "AGT", "TCA", "TCC", "TCG", "TCT"],
        "rna_codons": ["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"]
    },
    "T": {
        "single_letter": "T",
        "abbreviation": "Thr",
        "name": "Threonine",
        "dna_codons": ["ACA", "ACC", "ACG", "ACT"],
        "rna_codons": ["ACA", "ACC", "ACG", "ACU"]
    },
    "V": {
        "single_letter": "V",
        "abbreviation": "Val",
        "name": "Valine",
        "dna_codons": ["GTA", "GTC", "GTG", "GTT"],
        "rna_codons": ["GUA", "GUC", "GUG", "GUU"]
    },
    "W": {
        "single_letter": "W",
        "abbreviation": "Trp",
        "name": "Tryptophan",
        "dna_codons": ["TGG"],
        "rna_codons": ["UGG"]
    },
    "Y": {
        "single_letter": "Y",
        "abbreviation": "Tyr",
        "name": "Tyrosine",
        "dna_codons": ["TAC", "TAT"],
        "rna_codons": ["UAC", "UAU"]
    },
    "*": {
        "single_letter": "*",
        "abbreviation": "Stop",
        "name": "Stop Codon",
        "dna_codons": ["TAA", "TAG", "TGA"],
        "rna_codons": ["UAA", "UAG", "UGA"]
    }
}
