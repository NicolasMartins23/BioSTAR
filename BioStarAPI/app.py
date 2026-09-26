from __future__ import annotations

from fastapi import FastAPI, Request
from sqlalchemy import text

from BioStarAPI.auth.routes import router as auth_router
from BioStarAPI.auth.service import _check_short_rate_limit, enforce_request_limits, resolve_api_key
from BioStarAPI.controllers.analysis import (
    dna_to_protein,
    dna_to_rna,
    mutation_compare,
    protein_analysis,
    rna_to_dna,
    rna_to_protein,
)
from BioStarAPI.controllers.schemas import MutationCompareRequest, ProteinAnalysisRequest
from BioStarAPI.database.connection import engine

app = FastAPI(title="BioSTAR API", version="0.5.0")
app.state.rate_limit_state = {}


@app.middleware("http")
async def api_rate_limit_middleware(request: Request, call_next):
    path = request.url.path

    # Documentation and health checks are local operational endpoints and are not quota limited.
    if path.startswith("/api/auth/"):
        _check_short_rate_limit(request, authenticated=True)
    elif path.startswith("/api/"):
        api_key_id = resolve_api_key(request.headers.get("X-API-Key"))
        enforce_request_limits(request, api_key_id)

    return await call_next(request)


@app.get("/health")
def health() -> dict[str, str]:
    with engine.connect() as connection:
        connection.execute(text("SELECT 1"))
    return {"status": "ok"}


@app.get("/api/dna-rna")
def get_dna_rna(sequence: str) -> dict[str, str]:
    return dna_to_rna(sequence)


@app.get("/api/dna-protein")
def get_dna_protein(sequence: str) -> dict[str, str]:
    return dna_to_protein(sequence)


@app.get("/api/rna-protein")
def get_rna_protein(sequence: str) -> dict[str, str]:
    return rna_to_protein(sequence)


@app.get("/api/rna-dna")
def get_rna_dna(sequence: str) -> dict[str, str]:
    return rna_to_dna(sequence)


@app.post("/api/protein")
def post_protein(request: ProteinAnalysisRequest) -> dict[str, object]:
    return protein_analysis(request)


@app.post("/api/mutation_compare")
def post_mutation_compare(request: MutationCompareRequest) -> dict[str, object]:
    return mutation_compare(request.reference, request.sequence)


app.include_router(auth_router)
