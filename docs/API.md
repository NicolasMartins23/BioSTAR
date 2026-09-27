# BioSTAR API

BioSTAR provides targeted DNA, RNA and protein analysis through a FastAPI service.

## 1. Base URL

For a local Docker installation:

```text
http://localhost:8000
```

For a deployed installation, replace the host with the address provided by the BioSTAR administrator.

Interactive documentation:

```text
GET /docs
GET /redoc
```

Swagger UI at `/docs` is the recommended interactive reference because it is generated from the same OpenAPI schema used by the running API.

---

## 2. Authentication

BioSTAR supports anonymous requests with stricter rate limits and authenticated requests using an API key.

API keys are **not generated publicly**. Request an API key by emailing:

**api@bioapps.org**

Include:

- your name or organization;
- intended application or project;
- expected request volume;
- whether requests are interactive or automated.

After an API key is issued, send it in every authenticated request:

```http
X-API-Key: bst_live_...
```

Example:

```bash
curl \
  -H "X-API-Key: bst_live_YOUR_KEY" \
  "http://localhost:8000/api/dna-protein?sequence=ATGGCC"
```

### Security

Treat an API key as a password. Do not:

- commit it to Git;
- place it in a public JavaScript bundle;
- expose it in browser source code;
- put it in URLs;
- share it in issue trackers or logs.

For browser applications, requests requiring a private API key should normally go through your own backend rather than directly from the browser.

---

## 3. Rate limits

### Anonymous requests

Anonymous API requests are limited to:

- **1 request every 10 seconds per IP address**;
- **8,640 requests per server day in total**, regardless of the number of IP addresses.

### Authenticated requests

Authenticated API requests are limited to:

- **1 request every 2 seconds per IP address**;
- **86,400 requests per API key per server day**.

Authentication/key-management routes have their own short-window protection.

When a limit is exceeded, clients should expect HTTP `429 Too Many Requests` and back off before retrying.

Example client behavior:

```text
request
  ↓
200 OK → process response
  ↓
429 → wait → retry
```

Applications should avoid aggressive retry loops.

---

## 4. Sequence formats

BioSTAR normalizes sequence input by trimming whitespace and converting letters to uppercase.

DNA accepts:

```text
A C G T
```

RNA accepts:

```text
A C G U
```

Protein analysis accepts the standard amino-acid alphabet:

```text
A C D E F G H I K L M N P Q R S T V W Y
```

DNA/RNA conversion endpoints also accept a single FASTA record. Protein analysis requires exactly one FASTA record.

Example FASTA:

```fasta
>example_protein
MVLSPADKTNVKAAW
```

Multiple FASTA records are rejected by endpoints that require one sequence.

---

## 5. GET conversion endpoints

GET endpoints are designed for short, URL-friendly interactive conversions.

Maximum sequence length: **1,000 nucleotides**.

### 5.1 DNA → RNA

```http
GET /api/dna-rna?sequence=ATGGCC
```

cURL:

```bash
curl "http://localhost:8000/api/dna-rna?sequence=ATGGCC"
```

Response:

```json
{
  "sequence": "AUGGCC"
}
```

### 5.2 DNA → protein

```http
GET /api/dna-protein?sequence=ATGGCC
```

cURL:

```bash
curl "http://localhost:8000/api/dna-protein?sequence=ATGGCC"
```

Response:

```json
{
  "sequence": "MA"
}
```

### 5.3 RNA → protein

```http
GET /api/rna-protein?sequence=AUGGCC
```

Response:

```json
{
  "sequence": "MA"
}
```

### 5.4 RNA → DNA

```http
GET /api/rna-dna?sequence=AUGGCC
```

Response:

```json
{
  "sequence": "ATGGCC"
}
```

### GET error example

If a sequence exceeds 1,000 nucleotides:

```json
{
  "detail": "DNA sequence cannot exceed 1000 nucleotides"
}
```

If invalid characters are supplied:

```json
{
  "detail": "Invalid DNA sequence characters: X"
}
```

---

## 6. Protein analysis

Endpoint:

```http
POST /api/protein
Content-Type: application/json
```

The request contains one FASTA protein sequence and boolean flags selecting analyses.

### Request structure

```json
{
  "sequence": ">example\nMVLSPADKTNVKAAW",
  "get_full_test_results": false,
  "get_aminoacids_count": true,
  "get_isoelectric_point": true,
  "get_charge_at_pH": 7.0,
  "get_aromaticity": true,
  "get_secondary_structure_propensity": true,
  "get_molecular_weight": true,
  "get_hydrophobic_index": true,
  "get_composition_ratio": true,
  "get_extinction_coefficient": true
}
```

All analysis flags default to `false`.

### Full analysis

To request every available protein test:

```json
{
  "sequence": ">example\nMVLSPADKTNVKAAW",
  "get_full_test_results": true
}
```

When `get_full_test_results` is `true`, the individual test flags do not need to be supplied.

### Charge at pH

`get_charge_at_pH` is a strict floating-point value.

Valid:

```json
{
  "sequence": ">example\nMVLSPADKTNVKAAW",
  "get_charge_at_pH": 7.0
}
```

Invalid types such as a string should not be used:

```json
{
  "get_charge_at_pH": "7.0"
}
```

### cURL example

```bash
curl -X POST \
  "http://localhost:8000/api/protein" \
  -H "Content-Type: application/json" \
  -H "X-API-Key: bst_live_YOUR_KEY" \
  -d '{
    "sequence": ">example\nMVLSPADKTNVKAAW",
    "get_aminoacids_count": true,
    "get_isoelectric_point": true,
    "get_charge_at_pH": 7.0,
    "get_aromaticity": true
  }'
```

### Response structure

A response always includes the normalized sequence and its length:

```json
{
  "sequence": "MVLSPADKTNVKAAW",
  "length": 15,
  "aminoacids_count": {
    "M": 1,
    "V": 1,
    "L": 1
  },
  "isoelectric_point": 6.8,
  "charge_at_pH": {
    "pH": 7.0,
    "charge": -0.2
  },
  "aromaticity": 0.1333
}
```

Numeric values above are illustrative response-shape examples; the API calculates the actual values for the supplied sequence.

Possible analysis result fields are:

| Field | Meaning |
|---|---|
| `sequence` | Normalized protein sequence |
| `length` | Number of residues |
| `aminoacids_count` | Amino-acid counts |
| `isoelectric_point` | Estimated isoelectric point |
| `charge_at_pH` | Calculated charge and requested pH |
| `aromaticity` | Aromatic residue proportion |
| `secondary_structure_propensity` | Secondary-structure propensity data |
| `molecular_weight` | Molecular weight |
| `hydrophobic_index` | Hydrophobicity-related result |
| `composition_ratio` | Composition ratios |
| `extinction_coefficient` | Extinction coefficient result |

Only requested tests are included unless `get_full_test_results` is enabled.

---

## 7. Mutation comparison

Endpoint:

```http
POST /api/mutation_compare
Content-Type: application/json
```

This endpoint compares two coding DNA sequences.

### Request

```json
{
  "reference": "ATGGCCGAA",
  "sequence": "ATGGTCGAA"
}
```

cURL:

```bash
curl -X POST \
  "http://localhost:8000/api/mutation_compare" \
  -H "Content-Type: application/json" \
  -d '{
    "reference": "ATGGCCGAA",
    "sequence": "ATGGTCGAA"
  }'
```

### Requirements

Both sequences must:

1. contain valid DNA characters;
2. have the same length;
3. have a length divisible by three;
4. contain no more than 10,000 nucleotides.

### Response structure

```json
{
  "reference": "ATGGCCGAA",
  "sequence": "ATGGTCGAA",
  "mutations": []
}
```

The `mutations` field contains the mutation comparison produced by BioSTAR for the supplied sequences.

If lengths differ:

```json
{
  "detail": "Reference and sequence must have the same length"
}
```

If the length is not divisible by three:

```json
{
  "detail": "Reference and sequence lengths must be multiples of 3"
}
```

---

## 8. Batch conversion endpoints

Batch conversion routes use POST requests so that larger sequences do not need to be placed in a URL.

Maximum sequence length: **10,000 nucleotides**.

Available routes:

```text
POST /api/batch/dna-rna
POST /api/batch/dna-protein
POST /api/batch/rna-protein
POST /api/batch/rna-dna
```

### Request structure

All four endpoints use:

```json
{
  "sequence": "ATGGCC"
}
```

### Example

```bash
curl -X POST \
  "http://localhost:8000/api/batch/dna-protein" \
  -H "Content-Type: application/json" \
  -H "X-API-Key: bst_live_YOUR_KEY" \
  -d '{"sequence":"ATGGCC"}'
```

Response:

```json
{
  "sequence": "MA"
}
```

The batch endpoints currently represent the POST form of the sequence conversion operations; they are intended for larger payloads and future batch-oriented extensions.

---

## 9. HTTP status codes

Common responses include:

### `200 OK`

The operation completed successfully.

### `401 Unauthorized`

The request requires authentication or supplied credentials are not valid.

### `422 Unprocessable Entity`

The request structure or sequence content is invalid.

Examples include:

- empty sequence;
- invalid nucleotide/amino-acid character;
- multiple FASTA records where one is required;
- sequence exceeding the endpoint limit;
- mutation sequences with different lengths;
- mutation sequences whose lengths are not divisible by three.

### `429 Too Many Requests`

The applicable rate limit has been exceeded.

Clients should wait before retrying.

### `500 Internal Server Error`

An unexpected server-side error occurred. Contact the BioSTAR administrator if the problem persists.

---

## 10. Health check

The service exposes:

```http
GET /health
```

Example:

```bash
curl http://localhost:8000/health
```

Response:

```json
{
  "status": "ok"
}
```

The health endpoint verifies that the API can connect to PostgreSQL.

---

## 11. Recommended integration pattern

For a server-side application:

```text
Your application
      |
      | X-API-Key
      v
   BioSTAR API
      |
      v
  BioSTAR engine
```

For a browser application:

```text
Browser
   |
   v
Your backend
   |
   | X-API-Key
   v
BioSTAR API
```

Avoid putting a private BioSTAR API key directly into frontend JavaScript that is delivered to users.

---

## 12. JavaScript example

```javascript
const response = await fetch("http://localhost:8000/api/protein", {
  method: "POST",
  headers: {
    "Content-Type": "application/json",
    "X-API-Key": process.env.BIOSTAR_API_KEY,
  },
  body: JSON.stringify({
    sequence: ">example\nMVLSPADKTNVKAAW",
    get_full_test_results: true,
  }),
});

if (!response.ok) {
  throw new Error(`BioSTAR API error: ${response.status}`);
}

const result = await response.json();
console.log(result);
```

---

## 13. Python example

```python
import os
import requests

response = requests.post(
    "http://localhost:8000/api/protein",
    headers={"X-API-Key": os.environ["BIOSTAR_API_KEY"]},
    json={
        "sequence": ">example\nMVLSPADKTNVKAAW",
        "get_full_test_results": True,
    },
    timeout=30,
)
response.raise_for_status()
print(response.json())
```

---

## 14. OpenAPI / Swagger

For the authoritative machine-readable contract, use the running FastAPI service:

```text
/docs
/openapi.json
/redoc
```

The OpenAPI schema contains endpoint parameters, request models and interactive request controls.

When integrating BioSTAR into another application, prefer the OpenAPI schema over copying examples from this document because the running schema is the source of truth for the deployed API version.

---

## 15. API support

API key requests and API integration questions:

**api@bioapps.org**

When reporting an integration issue, include the endpoint, HTTP status code, sanitized request structure, and BioSTAR API version. Never include your API key in support requests.
