from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class AminoAcidData:
    symbol: str
    molecular_weight: float
    hydrophobicity: float
    alpha_helix: float
    beta_sheet: float
    pka: float | None
    pkb: float | None
    pkr: float | None


@dataclass(frozen=True)
class BiochemistryData:
    amino_acids: dict[str, AminoAcidData]
    aromatic: frozenset[str]
    nonpolar: frozenset[str]
    polar: frozenset[str]
    positive: frozenset[str]
    negative: frozenset[str]
    dna_codons: dict[str, str]
    rna_codons: dict[str, str]
    dna_stop_codons: frozenset[str]
    rna_stop_codons: frozenset[str]
    dna_start_codon: str
    rna_start_codon: str
    water_mass: float
    n_term_pka: float
    c_term_pka: float
