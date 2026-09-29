# BioSTAR

BioSTAR is an object-oriented Python library for biological sequence analysis.

## Architecture

The **BioSTAR package is the engine**. It contains the biological domain model,
reusable analysis functionality, and its immutable biochemical reference data.
The engine has no database dependency and does not use SQLite or any other
runtime database.

- `BioStar/NucleicAcids/` — DNA, RNA and nucleic-acid sequence functionality.
- `BioStar/Protein/` — protein sequence and biochemical analysis functionality.
- `BioStar/engine/reference_values.py` — bundled biochemical reference data as Python constants.
- `BioStar/analysis/` — sequence-analysis algorithms.
- `BioStar/io/` — input parsing such as FASTA.
- `BioStar/utils/` — framework-independent helpers.

The API is a separate application layer:

```
HTTP
  ↓
BioStarAPI/controllers
  ↓
BioStarAPI/services
  ↓
PostgreSQL
  ↓
API-owned application data

BioStarAPI/services
  ↓
BioSTAR engine
  ↓
Immutable Python reference data
```

## API foundation

`BioStarAPI/` contains the application layer for authentication, API keys,
request quotas, and HTTP endpoints.

PostgreSQL is used for API-owned application data. The Docker development and
production stacks start PostgreSQL as a dedicated service and run Alembic
migrations before starting the API.

Alembic migrations live in `BioStarAPI/database/migrations/`.

For a local installation, set `BIOSTAR_DATABASE_URL` to a PostgreSQL URL using
the `postgresql+psycopg://` SQLAlchemy driver.

## Public API

```python
from BioStar import DNA, RNA, Protein, NucleicAcid, OpenReadFrame
```

Biological features are organized by molecule type so that the engine remains
readable for biochemists and other life-science users.

## Docker

Development and production use PostgreSQL for API persistence. The BioSTAR
engine itself remains database-free; its biochemical reference data ships as
normal Python code.

```bash
docker compose -f docker-compose.dev.yml up --build
```

No biochemical seed step is required.
