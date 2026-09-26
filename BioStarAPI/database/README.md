# BioSTAR Database

The API owns database access. The BioSTAR engine does not import SQLAlchemy,
PostgreSQL, or FastAPI.

## Layers

```
Controller -> Service -> Repository -> SQLAlchemy -> PostgreSQL
                    |
                    -> BioSTAR engine
```

## Reference data

- amino_acids
- amino_acid_classes
- amino_acid_class_members
- nucleotides
- genetic_codes
- codons
- organisms
- codon_usage
- reference_sources

## Authentication

- users
- refresh_tokens

## Local database

Start PostgreSQL:

```bash
docker compose -f docker-compose.database.yml up -d
```

Set the connection string:

```bash
export BIOSTAR_DATABASE_URL="postgresql+psycopg://biostar:biostar@localhost:5432/biostar"
```

Create the schema:

```bash
alembic upgrade head
```

Load the initial biochemical reference data:

```bash
python -m BioStarAPI.database.seed
```

Schema migrations and scientific data imports remain separate.
