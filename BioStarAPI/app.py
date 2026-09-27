from __future__ import annotations

from fastapi import Depends, FastAPI, Request
from sqlalchemy import text

from BioStarAPI.auth.routes import router as auth_router
from BioStarAPI.auth.service import _check_short_rate_limit, enforce_request_limits, resolve_api_key
from BioStarAPI.controllers.analysis import GET_SEQUENCE_MAX_LENGTH, POST_SEQUENCE_MAX_LENGTH
from BioStarAPI.controllers.schemas import BatchSequenceRequest, MutationCompareRequest, ProteinAnalysisRequest
from BioStarAPI.database.connection import engine
from BioStarAPI.database.repositories.biochemistry import BiochemistryRepository
from BioStarAPI.openapi import custom_openapi
from BioStarAPI.resources.exceptions import BioStarAPIError, api_error_response, register_exception_handlers
from BioStarAPI.resources.responses import APIResponse
from BioStarAPI.services.dependencies import get_biochemistry_repository
from BioStarAPI.services.mutation_service import MutationService
from BioStarAPI.services.protein_service import ProteinService
from BioStarAPI.services.sequence_service import SequenceService

app = FastAPI(title="BioSTAR API", version="0.5.0")
app.state.rate_limit_state = {}
register_exception_handlers(app)


def get_sequence_service(
    repository: BiochemistryRepository = Depends(get_biochemistry_repository),
) -> SequenceService:
    return SequenceService(repository)


@app.middleware("http")
async def api_rate_limit_middleware(request: Request, call_next):
    path = request.url.path

    if path == "/health" or path == "/api/v1/health":
        return await call_next(request)

    try:
        if path.startswith("/api/auth/"):
            _check_short_rate_limit(request, authenticated=True)
        elif path.startswith("/api/"):
            api_key_id = resolve_api_key(request.headers.get("X-API-Key"))
            enforce_request_limits(request, api_key_id)

        return await call_next(request)
    except BioStarAPIError as exception:
        return api_error_response(exception)


def _check_database_health() -> dict[str, str]:
    with engine.connect() as connection:
        connection.execute(text("SELECT 1"))
    return {"status": "ok"}


@app.get("/health", tags=["System"], summary="Check API health")
def health() -> dict[str, str]:
    return _check_database_health()


@app.get("/api/v1/health", tags=["System"], summary="Check API health")
def api_health() -> APIResponse[dict[str, str]]:
    return APIResponse(data=_check_database_health())


@app.get("/api/dna-rna", tags=["Conversions"])
def get_dna_rna(
    sequence: str,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.dna_to_rna(sequence, GET_SEQUENCE_MAX_LENGTH))


@app.get("/api/dna-protein", tags=["Conversions"])
def get_dna_protein(
    sequence: str,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.dna_to_protein(sequence, GET_SEQUENCE_MAX_LENGTH))


@app.get("/api/rna-protein", tags=["Conversions"])
def get_rna_protein(
    sequence: str,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.rna_to_protein(sequence, GET_SEQUENCE_MAX_LENGTH))


@app.get("/api/rna-dna", tags=["Conversions"])
def get_rna_dna(
    sequence: str,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.rna_to_dna(sequence, GET_SEQUENCE_MAX_LENGTH))


@app.post("/api/protein", tags=["Protein"])
def post_protein(
    request: ProteinAnalysisRequest,
    repository: BiochemistryRepository = Depends(get_biochemistry_repository),
) -> APIResponse[dict[str, object]]:
    return APIResponse(data=ProteinService(repository).analyze(request))


@app.post("/api/mutation_compare", tags=["Mutations"])
def post_mutation_compare(
    request: MutationCompareRequest,
    repository: BiochemistryRepository = Depends(get_biochemistry_repository),
) -> APIResponse[dict[str, object]]:
    return APIResponse(data=MutationService(repository).compare(request.reference, request.sequence))


@app.post("/api/batch/dna-rna", tags=["Batch"])
def batch_dna_rna(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.dna_to_rna(request.sequence, POST_SEQUENCE_MAX_LENGTH))


@app.post("/api/batch/dna-protein", tags=["Batch"])
def batch_dna_protein(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.dna_to_protein(request.sequence, POST_SEQUENCE_MAX_LENGTH))


@app.post("/api/batch/rna-protein", tags=["Batch"])
def batch_rna_protein(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]:
    return APIResponse(data=service.rna_to_protein(request.sequence, POST_SEQUENCE_MAX_LENGTH))


@app.post("/api/batch/rna-dna", tags=["Batch"])
def batch_rna_dna(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> APIResponse[dict[str, str]]
    return APIResponse(data=service.rna_to_dna(request.sequence, POST_SEQUENCE_MAX_LENGTH))


app.include_router(auth_router)
app.openapi = lambda: custom_openapi(app)
