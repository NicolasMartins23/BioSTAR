from __future__ import annotations

from fastapi import FastAPI, Request
from sqlalchemy import text

from BioStarAPI.auth.routes import router as auth_router
from BioStarAPI.auth.service import _check_short_rate_limit, enforce_request_limits, resolve_api_key
from BioStarAPI.controllers.analysis import (
    batch_dna_protein,
    batch_dna_rna,
    batch_rna_dna,
    batch_rna_protein,
    dna_to_protein,
    dna_to_rna,
    mutation_compare,
    protein_analysis,
    rna_to_dna,
    rna_to_protein,
)
from BioStarAPI.database.connection import engine
from BioStarAPI.openapi import custom_openapi
from BioStarAPI.resources.exceptions import BioStarAPIError, api_error_response, register_exception_handlers
from BioStarAPI.resources.messages import MessageCode
from BioStarAPI.resources.responses import APIResponse

app = FastAPI(title="BioSTAR API", version="0.5.0")
app.state.rate_limit_state = {}
register_exception_handlers(app)


@app.middleware("http")
async def api_rate_limit_middleware(request: Request, call_next):
    path = request.url.path

    if path == "/health" or path == "/api/v1/health":
        return await call_next(request)

    try:
        if path.startswith("/api/auth/"):
            _check_short_rate_limit(request, authenticated=True)
        elif path.startswith("/api/"):
            raw_api_key = request.headers.get("X-API-Key")
            api_key_id = resolve_api_key(raw_api_key)
            if raw_api_key is not None and api_key_id is None:
                raise BioStarAPIError(401, MessageCode.INVALID_API_KEY)
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


app.add_api_route("/api/dna-rna", dna_to_rna, methods=["GET"], tags=["Conversions"])
app.add_api_route("/api/dna-protein", dna_to_protein, methods=["GET"], tags=["Conversions"])
app.add_api_route("/api/rna-protein", rna_to_protein, methods=["GET"], tags=["Conversions"])
app.add_api_route("/api/rna-dna", rna_to_dna, methods=["GET"], tags=["Conversions"])
app.add_api_route("/api/protein", protein_analysis, methods=["POST"], tags=["Protein"])
app.add_api_route("/api/mutation_compare", mutation_compare, methods=["POST"], tags=["Mutations"])
app.add_api_route("/api/batch/dna-rna", batch_dna_rna, methods=["POST"], tags=["Batch"])
app.add_api_route("/api/batch/dna-protein", batch_dna_protein, methods=["POST"], tags=["Batch"])
app.add_api_route("/api/batch/rna-protein", batch_rna_protein, methods=["POST"], tags=["Batch"])
app.add_api_route("/api/batch/rna-dna", batch_rna_dna, methods=["POST"], tags=["Batch"])

app.include_router(auth_router)
app.openapi = lambda: custom_openapi(app)
