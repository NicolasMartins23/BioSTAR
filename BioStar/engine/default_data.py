from __future__ import annotations

from functools import lru_cache

from BioStar.engine.biochemistry import BiochemistryData
from BioStar.engine.reference_values import (
    AMINO_ACIDS,
    AROMATIC,
    C_TERM_PKA,
    DNA_CODONS,
    DNA_START_CODON,
    DNA_STOP_CODONS,
    NEGATIVE,
    NONPOLAR,
    N_TERM_PKA,
    POLAR,
    POSITIVE,
    RNA_CODONS,
    RNA_START_CODON,
    RNA_STOP_CODONS,
    WATER_MASS,
)


@lru_cache(maxsize=1)
def get_default_biochemistry() -> BiochemistryData:
    return BiochemistryData(
        amino_acids=dict(AMINO_ACIDS),
        aromatic=AROMATIC,
        nonpolar=NONPOLAR,
        polar=POLAR,
        positive=POSITIVE,
        negative=NEGATIVE,
        dna_codons=dict(DNA_CODONS),
        rna_codons=dict(RNA_CODONS),
        dna_stop_codons=DNA_STOP_CODONS,
        rna_stop_codons=RNA_STOP_CODONS,
        dna_start_codon=DNA_START_CODON,
        rna_start_codon=RNA_START_CODON,
        water_mass=WATER_MASS,
        n_term_pka=N_TERM_PKA,
        c_term_pka=C_TERM_PKA,
    )
