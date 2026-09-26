from __future__ import annotations

from fastapi import FastAPI
from sqlalchemy import text

from BioStarAPI.controllers.analysis import (
    GET_SEQUENCE_MAX_LENGTH,
    dna_to_protein,
    dna_to_rna,
    mutation_compare,
    protein_analysis,
    rna_to_dna,
    rna_to_protein,
)
from BioStarAPI.controllers.schemas import BatchSequenceRequest, MutationCompareRequest, ProteinAnalysisRequest
from BioStarAPI.database.connection import engine

app = FastAPI(title="BioSTAR API", version="0.4.0")


@app.get("/health")
def health() -> dict[str, str]:
    with engine.connect() as connection:
        connection.execute(text("SELECT 1"))
    return {"status": "ok"}


@app.get("/api/dna-rna")
def get_dna_rna(sequence: str) -> dict[str, str]:
    return dna_to_rna(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.get("/api/dna-protein")
def get_dna_protein(sequence: str) -> dict[str, str]:
    return dna_to_protein(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.get("/api/rna-protein")
def get_rna_protein(sequence: str) -> dict[str, str]:
    return rna_to_protein(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.get("/api/rna-dna")
def get_rna_dna(sequence: str) -> dict[str, str]:
    return rna_to_dna(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.post("/api/protein")
def post_protein(request: ProteinAnalysisRequest) -> dict[str, object]:
    return protein_analysis(request)


@app.post("/api/mutation_compare")
def post_mutation_compare(request: MutationCompareRequest) -> dict[str, object]:
    return mutation_compare(request.reference, request.sequence)


@app.post("/api/batch/dna-rna")
def batch_dna_rna(request: BatchSequenceRequest) -> dict[str, str]:
    return dna_to_rna(request.sequence)


@app.post("/api/batch/dna-protein")
def batch_dna_protein(request: BatchSequenceRequest) -> dict[str, str]:
    return dna_to_protein(request.sequence)


@app.post("/api/batch/rna-protein")
def batch_rna_protein(request: BatchSequenceRequest) -> dict[str, str]:
    return rna_to_protein(request.sequence)


@app.post("/api/batch/rna-dna")
def batch_rna_dna(request: BatchSequenceRequest) -> dict[str, str]:
    return rna_to_dna(request.sequence)
