from fastapi.openapi.utils import get_openapi


def custom_openapi(app):
    if app.openapi_schema:
        return app.openapi_schema

    schema = get_openapi(
        title="BioSTAR API",
        version="0.5.0",
        description="""
BioSTAR is a bioinformatics API for targeted DNA, RNA and protein analysis.

## Authentication

Interactive and public requests may be made without an API key subject to anonymous rate limits.
For higher-volume programmatic access, request an API key by email at **api@bioapps.org**.

Send requests to `api@bioapps.org` with a short description of your intended use and expected request volume.

Authenticated requests use:

`X-API-Key: bst_live_...`

## Sequence limits

- GET conversion endpoints: maximum 1,000 nucleotides.
- POST analysis endpoints: maximum 10,000 sequence characters.
- Batch endpoints: maximum 10,000 nucleotides per sequence.

## Formats

Protein analysis accepts a single FASTA sequence. DNA/RNA conversion endpoints accept raw sequences.
""",
        routes=app.routes,
    )

    schema["info"]["contact"] = {
        "name": "BioSTAR API Support",
        "email": "api@bioapps.org",
    }
    schema["servers"] = [{"url": "/", "description": "Current BioSTAR server"}]
    schema["tags"] = [
        {"name": "Conversions", "description": "DNA/RNA/protein sequence conversions."},
        {"name": "Protein", "description": "Protein sequence analysis."},
        {"name": "Mutations", "description": "Mutation comparison."},
        {"name": "Batch", "description": "POST endpoints for larger sequence requests."},
        {"name": "Authentication", "description": "API key access and authentication."},
        {"name": "System", "description": "Service health and operational endpoints."},
    ]

    app.openapi_schema = schema
    return schema
