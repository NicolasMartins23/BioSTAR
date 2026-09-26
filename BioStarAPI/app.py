from __future__ import annotations

from fastapi import FastAPI, Request
from sqlalchemy import text

from BioStarAPI.auth.routes import router as auth_router
from BioStarAPI.auth.service import _check_short_rate_limit, enforce_request_limits, resolve_api_key
from BioStarAPI.controllers.analysis import (
    GET_SEQUENCE_MAX_LENGTH,
    POST_SEQUENCE_MAX_LENGTH,
)
from BioStarAPI.controllers.schemas import BatchSequenceRequest, MutationCompareRequest, ProteinAnalysisRequest
from BioStarAPI.database.connection import engine
from BioStarAPI.database.repositories.biochemistry import BiochemistryRepository
from BioStarAPI.services.dependencies import get_biochemistry_repository
from BioStarAPI.services.mutation_service import MutationService
from BioStarAPI.services.protein_service import ProteinService
from BioStarAPI.services.sequence_service import SequenceService
from fastapi import Depends
from BioStarAPI.openapi import custom_openapi

app = FastAPI(title="BioSTAR API", version="0.5.0")
app.state.rate_limit_state = {}


def _sequence_example(length: int = 24) -> str:
    return ("ATGGCCGAACTGCTGATCGTTAC" * ((length + 23) // 24))[:length]


@app.middleware("http")
async def api_rate_limit_middleware(request: Request, call_next):
    path = request.url.path
    if path.startswith("/api/auth/"):
        _check_short_rate_limit(request, authenticated=True)
    elif path.startswith("/api/"):
        api_key_id = resolve_api_key(request.headers.get("X-API-Key"))
        enforce_request_limits(request, api_key_id)
    return await call_next(request)


@app.get("/health", tags=["System"], summary="Check API health")
def health() -> dict[str, str]:
    """Return a simple health status after verifying the PostgreSQL connection."""
    with engine.connect() as connection:
        connection.execute(text("SELECT 1"))
    return {"status": "ok"}


@app.get(
    "/api/dna-rna",
    tags=["Conversions"],
    summary="Transcribe DNA to RNA",
    description=f"Convert a DNA sequence to RNA. Intended for short interactive requests. Maximum {GET_SEQUENCE_MAX_LENGTH:,} nucleotides.",
    response_description="The resulting RNA sequence.",
)
def get_dna_rna(
    sequence: str,
) -> dict[str, str]:
    """Convert DNA bases to their RNA representation."""
    return dna_to_rna(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.get(
    "/api/dna-protein",
    tags=["Conversions"],
    summary="Translate DNA to protein",
    description=f"Translate a DNA sequence into its protein sequence. Maximum {GET_SEQUENCE_MAX_LENGTH:,} nucleotides.",
    response_description="The resulting protein sequence.",
)
def get_dna_protein(sequence: str) -> dict[str, str]:
    """Translate DNA into a protein sequence."""
    return dna_to_protein(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.get(
    "/api/rna-protein",
    tags=["Conversions"],
    summary="Translate RNA to protein",
    description=f"Translate an RNA sequence into its protein sequence. Maximum {GET_SEQUENCE_MAX_LENGTH:,} nucleotides.",
    response_description="The resulting protein sequence.",
)
def get_rna_protein(sequence: str) -> dict[str, str]:
    """Translate RNA into a protein sequence."""
    return rna_to_protein(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.get(
    "/api/rna-dna",
    tags=["Conversions"],
    summary="Convert RNA to DNA",
    description=f"Convert an RNA sequence to DNA. Maximum {GET_SEQUENCE_MAX_LENGTH:,} nucleotides.",
    response_description="The resulting DNA sequence.",
)
def get_rna_dna(sequence: str) -> dict[str, str]:
    """Convert RNA bases to their DNA representation."""
    return rna_to_dna(sequence, GET_SEQUENCE_MAX_LENGTH)


@app.post(
    "/api/protein",
    tags=["Protein"],
    summary="Analyze a protein sequence",
    description=(
        f"Run selected protein analyses on exactly one FASTA sequence. "
        f"The normalized protein sequence is limited to {POST_SEQUENCE_MAX_LENGTH:,} residues. "
        "Set get_full_test_results to true to request every available analysis."
    ),
    response_description="The requested protein analysis results.",
)
def post_protein(request: ProteinAnalysisRequest) -> dict[str, object]:
    """Run one or more supported protein analyses."""
    return protein_analysis(request)


@app.post(
    "/api/mutation_compare",
    tags=["Mutations"],
    summary="Compare two DNA sequences",
    description=(
        f"Compare a reference DNA sequence with a mutated sequence. Both sequences must have equal lengths, "
        f"be divisible by three, and contain at most {POST_SEQUENCE_MAX_LENGTH:,} nucleotides."
    ),
    response_description="The two normalized sequences and detected mutations.",
)
def post_mutation_compare(request: MutationCompareRequest, repository: BiochemistryRepository = Depends(get_biochemistry_repository)) -> dict[str, object]:
    """Compare two coding DNA sequences and return detected mutations."""
    return MutationService(repository).compare(request.reference, request.sequence)


@app.post(
    "/api/batch/dna-rna",
    tags=["Batch"],
    summary="Convert a larger DNA sequence to RNA",
    description=f"POST variant of DNA-to-RNA conversion. Maximum {POST_SEQUENCE_MAX_LENGTH:,} nucleotides.",
)
def batch_dna_rna(request: BatchSequenceRequest, service: SequenceService = Depends(get_sequence_service)) -> dict[str, str]:
    return service.dna_to_rna(request.sequence, POST_SEQUENCE_MAX_LENGTH)


@app.post(
    "/api/batch/dna-protein",
    tags=["Batch"],
    summary="Translate a larger DNA sequence to protein",
    description=f"POST variant of DNA-to-protein translation. Maximum {POST_SEQUENCE_MAX_LENGTH:,} nucleotides.",
)
def batch_dna_protein(request: BatchSequenceRequest, service: SequenceService = Depends(get_sequence_service)) -> dict[str, str]:
    return service.dna_to_protein(request.sequence, POST_SEQUENCE_MAX_LENGTH)


@app.post(
    "/api/batch/rna-protein",
    tags=["Batch"],
    summary="Translate a larger RNA sequence to protein",
    description=f"POST variant of RNA-to-protein translation. Maximum {POST_SEQUENCE_MAX_LENGTH:,} nucleotides.",
)
def batch_rna_protein(request: BatchSequenceRequest, service: SequenceService = Depends(get_sequence_service)) -> dict[str, str]:
    return service.rna_to_protein(request.sequence, POST_SEQUENCE_MAX_LENGTH)


@app.post(
    "/api/batch/rna-dna",
    tags=["Batch"],
    summary="Convert a larger RNA sequence to DNA",
    description=f"POST variant of RNA-to-DNA conversion. Maximum {POST_SEQUENCE_MAX_LENGTH:,} nucleotides.",
)
def batch_rna_dna(request: BatchSequenceRequest, service: SequenceService = Depends(get_sequence_service)) -> dict[str, str]:
    return service.rna_to_dna(request.sequence, POST_SEQUENCE_MAX_LENGTH)


app.include_router(auth_router)
app.openapi = lambda: custom_openapi(app)
