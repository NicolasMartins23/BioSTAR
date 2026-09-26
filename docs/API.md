# BioSTAR API

BioSTAR provides targeted DNA, RNA and protein analysis through a FastAPI service.

## Interactive documentation

When the server is running, open `/docs` for Swagger UI and `/redoc` for ReDoc.

## Authentication

API keys are not generated anonymously. Request a key by emailing **api@bioapps.org**. Include:

- intended application or project;
- expected request volume;
- whether requests are interactive or automated.

After approval, send the key using the `X-API-Key` header:

```http
X-API-Key: bst_live_...
```

Never commit an API key to source control or expose it in client-side code.

## Endpoints

### Conversions

```text
GET /api/dna-rna
GET /api/dna-protein
GET /api/rna-protein
GET /api/rna-dna
```

GET conversion requests are intended for short interactive sequences and accept up to 1,000 nucleotides.

### Protein analysis

```text
POST /api/protein
```

Accepts a FASTA protein sequence and optional analysis flags.

### Mutation comparison

```text
POST /api/mutation_compare
```

```json
{
  "reference": "ATGGCCGAA",
  "sequence": "ATGGTCGAA"
}
```

### Batch conversions

```text
POST /api/batch/dna-rna
POST /api/batch/dna-protein
POST /api/batch/rna-protein
POST /api/batch/rna-dna
```

Batch requests use POST bodies and support sequences up to 10,000 nucleotides.

## Rate limits

Unauthenticated requests are limited to one request per 10 seconds per IP and 8,640 requests per server day.

Authenticated requests are limited to one request per 2 seconds per IP and 86,400 requests per API key per server day.

Clients should treat HTTP `429` as a signal to back off and retry later.

## Example

```bash
curl -H "X-API-Key: bst_live_..." \
  "https://your-biostar-host/api/dna-protein?sequence=ATGGCC"
```

For the authoritative endpoint schema, parameters, request bodies and response models, use the generated Swagger/OpenAPI documentation at `/docs`.
