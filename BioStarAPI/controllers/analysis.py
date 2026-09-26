from __future__ import annotations

import re
from typing import Any

from BioStar.analysis.sequence_comparison import CompareNucleotideSequence
from BioStar.domain.nucleic_acid import DNA, RNA
from BioStar.domain.protein import Protein
from BioStar.io.fasta import FastaParser
from fastapi import HTTPException

GET_SEQUENCE_MAX_LENGTH = 1_000
POST_SEQUENCE_MAX_LENGTH = 10_000
DNA_ALPHABET = set("ACGT")
RNA_ALPHABET = set("ACGU")
PROTEIN_ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")


def _normalize_single_sequence(
    sequence: str,
    alphabet: set[str],
    kind: str,
    max_length: int,
) -> str:
    value = sequence.strip().upper()
    if value.startswith(">"):
        records = FastaParser().get_sequence_map(value)
        if len(records) != 1:
            raise HTTPException(status_code=422, detail="Exactly one FASTA sequence is required")
        value = str(records[0]["sequence"]).replace(" ", "").replace("\r", "").replace("\n", "")
    else:
        value = re.sub(r"\s+", "", value)

    if not value:
        raise HTTPException(status_code=422, detail=f"{kind} sequence cannot be empty")
    if len(value) > max_length:
        raise HTTPException(status_code=422, detail=f"{kind} sequence cannot exceed {max_length} nucleotides")
    invalid = sorted(set(value) - alphabet)
    if invalid:
        raise HTTPException(status_code=422, detail=f"Invalid {kind} sequence characters: {', '.join(invalid)}")
    return value


def normalize_dna(sequence: str, max_length: int = POST_SEQUENCE_MAX_LENGTH) -> str:
    return _normalize_single_sequence(sequence, DNA_ALPHABET, "DNA", max_length)


def normalize_rna(sequence: str, max_length: int = POST_SEQUENCE_MAX_LENGTH) -> str:
    return _normalize_single_sequence(sequence, RNA_ALPHABET, "RNA", max_length)


def normalize_protein_fasta(sequence: str) -> str:
    value = sequence.strip()
    records = FastaParser().get_sequence_map(value)
    if len(records) != 1:
        raise HTTPException(status_code=422, detail="Exactly one FASTA protein sequence is required")
    protein = re.sub(r"\s+", "", str(records[0]["sequence"]).upper())
    if not protein:
        raise HTTPException(status_code=422, detail="Protein sequence cannot be empty")
    if len(protein) > POST_SEQUENCE_MAX_LENGTH:
        raise HTTPException(status_code=422, detail=f"Protein sequence cannot exceed {POST_SEQUENCE_MAX_LENGTH} residues")
    invalid = sorted(set(protein) - PROTEIN_ALPHABET)
    if invalid:
        raise HTTPException(status_code=422, detail=f"Invalid protein sequence characters: {', '.join(invalid)}")
    return protein


def dna_to_rna(sequence: str, max_length: int = POST_SEQUENCE_MAX_LENGTH) -> dict[str, Any]:
    normalized = normalize_dna(sequence, max_length)
    return {"sequence": DNA(normalized).rna_sequence()}


def dna_to_protein(sequence: str, max_length: int = POST_SEQUENCE_MAX_LENGTH) -> dict[str, Any]:
    normalized = normalize_dna(sequence, max_length)
    protein = DNA(normalized).to_protein()
    return {"sequence": protein.sequence}


def rna_to_protein(sequence: str, max_length: int = POST_SEQUENCE_MAX_LENGTH) -> dict[str, Any]:
    normalized = normalize_rna(sequence, max_length)
    protein = RNA(normalized).to_protein()
    return {"sequence": protein.sequence}


def rna_to_dna(sequence: str, max_length: int = POST_SEQUENCE_MAX_LENGTH) -> dict[str, Any]:
    normalized = normalize_rna(sequence, max_length)
    return {"sequence": RNA(normalized).dna_sequence()}


def protein_analysis(request: Any) -> dict[str, Any]:
    sequence = normalize_protein_fasta(request.sequence)
    protein = Protein(sequence)

    requested = {
        "aminoacids_count": request.get_aminoacids_count,
        "isoelectric_point": request.get_isoelectric_point,
        "charge_at_pH": request.get_charge_at_pH is not None,
        "aromaticity": request.get_aromaticity,
        "secondary_structure_propensity": request.get_secondary_structure_propensity,
        "molecular_weight": request.get_molecular_weight,
        "hydrophobic_index": request.get_hydrophobic_index,
        "composition_ratio": request.get_composition_ratio,
        "extinction_coefficient": request.get_extinction_coefficient,
    }
    if request.get_full_test_results:
        requested = {name: True for name in requested}

    results: dict[str, Any] = {"sequence": protein.sequence, "length": protein.sequence_size}
    if requested["aminoacids_count"]:
        results["aminoacids_count"] = protein.count["by_aminoacid"]
    if requested["isoelectric_point"]:
        results["isoelectric_point"] = protein.isoelectric_point()
    if requested["charge_at_pH"]:
        pH = request.get_charge_at_pH if request.get_charge_at_pH is not None else 7.0
        results["charge_at_pH"] = {"pH": pH, "charge": protein.charge_at_pH(pH)}
    if requested["aromaticity"]:
        results["aromaticity"] = protein.aromacity()
    if requested["secondary_structure_propensity"]:
        results["secondary_structure_propensity"] = protein.secondary_structure_propensity()
    if requested["molecular_weight"]:
        results["molecular_weight"] = protein.molecular_weight()
    if requested["hydrophobic_index"]:
        results["hydrophobic_index"] = protein.hydrophobic_index()
    if requested["composition_ratio"]:
        results["composition_ratio"] = protein.composition_ratio()
    if requested["extinction_coefficient"]:
        results["extinction_coefficient"] = protein.extinction_coefficient()
    return results


def mutation_compare(reference: str, sequence: str) -> dict[str, Any]:
    reference_dna = normalize_dna(reference, POST_SEQUENCE_MAX_LENGTH)
    sequence_dna = normalize_dna(sequence, POST_SEQUENCE_MAX_LENGTH)
    if len(reference_dna) != len(sequence_dna):
        raise HTTPException(status_code=422, detail="Reference and sequence must have the same length")
    if len(reference_dna) % 3 != 0:
        raise HTTPException(status_code=422, detail="Reference and sequence lengths must be multiples of 3")
    mutations = CompareNucleotideSequence(reference_dna, sequence_dna).compare(show_only_mutations=True)
    return {"reference": reference_dna, "sequence": sequence_dna, "mutations": mutations}
