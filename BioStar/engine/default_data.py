from __future__ import annotations

from BioStar.data.biochemistry import (
    AMINOACID_TABLE,
    AMINOACIDS_AROMATIC,
    AMINOACIDS_NEGATIVE,
    AMINOACIDS_NONPOLAR,
    AMINOACIDS_POLAR,
    AMINOACIDS_POSITIVE,
    C_TERM_PKA,
    N_TERM_PKA,
    START_CODON_DNA,
    START_CODON_RNA,
    STOP_CODON_DNA,
    STOP_CODON_RNA,
    TABLE_DNA_CODON_TO_AMINOACID,
    TABLE_RNA_CODON_TO_AMINOACID,
    WATER_MASS,
)
from BioStar.engine.biochemistry import AminoAcidData, BiochemistryData


def get_default_biochemistry() -> BiochemistryData:
    return BiochemistryData(
        amino_acids={
            symbol: AminoAcidData(
                symbol=symbol,
                molecular_weight=float(data["weight"]),
                hydrophobicity=float(data["hydrophobicity"]),
                alpha_helix=float(data["alpha_helix"]),
                beta_sheet=float(data["beta_sheet"]),
                pka=float(data["pKa"]) if data.get("pKa") is not None else None,
                pkb=float(data["pKb"]) if data.get("pKb") is not None else None,
                pkr=float(data["pKr"]) if data.get("pKr") is not None else None,
            )
            for symbol, data in AMINOACID_TABLE.items()
            if symbol != "*"
        },
        aromatic=frozenset(AMINOACIDS_AROMATIC),
        nonpolar=frozenset(AMINOACIDS_NONPOLAR),
        polar=frozenset(AMINOACIDS_POLAR),
        positive=frozenset(AMINOACIDS_POSITIVE),
        negative=frozenset(AMINOACIDS_NEGATIVE),
        dna_codons=dict(TABLE_DNA_CODON_TO_AMINOACID),
        rna_codons=dict(TABLE_RNA_CODON_TO_AMINOACID),
        dna_stop_codons=frozenset(STOP_CODON_DNA),
        rna_stop_codons=frozenset(STOP_CODON_RNA),
        dna_start_codon=START_CODON_DNA,
        rna_start_codon=START_CODON_RNA,
        water_mass=WATER_MASS,
        n_term_pka=N_TERM_PKA,
        c_term_pka=C_TERM_PKA,
    )
