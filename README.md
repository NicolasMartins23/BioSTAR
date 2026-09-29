# BioSTAR

BioSTAR is an object-oriented Python library for biological sequence analysis.

## Architecture

The **BioSTAR package is the engine**. It contains the biological domain model,
reusable analysis functionality, and its own immutable biochemical reference
database. It does not depend on FastAPI, SQLAlchemy, authentication, or HTTP.

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
SQLite
  ↓
API-owned application data

BioStarAPI/services
  ↓
BioSTAR engine
  ↓
Bundled SQLite biochemical reference data
```

## API foundation

`BioStarAPI/` contains the application layer for authentication, API keys,
request quotas, and HTTP endpoints.

SQLite is the only database technology used by the project. API-owned
application data is stored in `biostar.db`; the Docker development and
production stacks persist it in a dedicated volume.

Alembic migrations live in `BioStarAPI/database/migrations/`.

For a local installation:

```bash
export BIOSTAR_DATABASE_URL="sqlite:///./biostar.db"
alembic upgrade head
```

## Public API

```python
from BioStar import DNA, RNA, Protein, NucleicAcid, OpenReadFrame
```

Biological features are organized by molecule type so that the engine remains
readable for biochemists and other life-science users.

## Docker

Development and production use the same SQLite database architecture. No
external database server or database administration service is required.

```bash
docker compose -f docker-compose.dev.yml up --build
```

No biochemical seed step is required.
