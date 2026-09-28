# BioSTAR

BioSTAR is an object-oriented Python library for biological sequence analysis.

## Architecture

The **BioSTAR package is the engine**. It contains the biological domain model,
reusable analysis functionality, and its own immutable biochemical reference
database. It does not depend on FastAPI, SQLAlchemy, PostgreSQL,
authentication, or HTTP.

- `BioStar/NucleicAcids/` — DNA, RNA and nucleic-acid sequence functionality.
- `BioStar/Protein/` — protein sequence and biochemical analysis functionality.
- `BioStar/data/` — bundled SQLite biochemical reference database.
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
SQLite biochemical reference data
```

## API foundation

`BioStarAPI/` contains the application layer for authentication, API keys,
request quotas, and HTTP endpoints.

PostgreSQL stores API-owned application state. Scientific reference data is
bundled with the BioSTAR engine instead of being seeded into PostgreSQL.

Alembic migrations live in `BioStarAPI/database/migrations/`.

Database migrations use:

```bash
export BIOSTAR_DATABASE_URL="postgresql+psycopg://user:password@localhost:5432/biostar"
alembic upgrade head
```

## Public API

```python
from BioStar import DNA, RNA, Protein, NucleicAcid, OpenReadFrame
```

Biological features are organized by molecule type so that the engine remains
readable for biochemists and other life-science users.

## Local PostgreSQL

The repository includes PostgreSQL for API development:

```bash
docker compose -f docker-compose.database.yml up -d
export BIOSTAR_DATABASE_URL="postgresql+psycopg://biostar:biostar@localhost:5432/biostar"
alembic upgrade head
```

No biochemical seed step is required.
