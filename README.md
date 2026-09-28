# BioSTAR

BioSTAR is an object-oriented Python library for biological sequence analysis.

## Architecture

The **BioSTAR package is the engine**. It contains the biological domain model and
reusable analysis functionality and does not depend on FastAPI, SQLAlchemy,
PostgreSQL, authentication, or HTTP.

- `BioStar/NucleicAcids/` — DNA, RNA and nucleic-acid sequence functionality.
- `BioStar/Protein/` — protein sequence and biochemical analysis functionality.
- `BioStar/data/` — current in-memory biochemical reference data.
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
BioStarAPI/database/repositories
  ↓
SQLAlchemy
  ↓
PostgreSQL

BioStarAPI/services
  ↓
BioSTAR engine
```

## API foundation

`BioStarAPI/` contains the initial database and application-layer structure.

The database currently models:

- amino acids and their biochemical properties
- amino acid classifications
- nucleotides
- genetic codes and codons
- organisms and codon usage
- scientific reference sources
- users and refresh tokens

Alembic migrations live in `BioStarAPI/database/migrations/`.

Database migrations use:

```bash
export BIOSTAR_DATABASE_URL="postgresql+psycopg://user:password@localhost:5432/biostar"
alembic upgrade head
```

Scientific reference data is intentionally not inserted by the initial schema
migration. It will be seeded separately so schema migrations and scientific
data imports remain independent.

## Public API

```python
from BioStar import DNA, RNA, Protein, NucleicAcid, OpenReadFrame
```

Biological features are organized by molecule type so that the engine remains
readable for biochemists and other life-science users.


## Local PostgreSQL

The repository includes a PostgreSQL development database:

```bash
docker compose -f docker-compose.database.yml up -d
export BIOSTAR_DATABASE_URL="postgresql+psycopg://biostar:biostar@localhost:5432/biostar"
alembic upgrade head
python -m BioStarAPI.database.seed
```
