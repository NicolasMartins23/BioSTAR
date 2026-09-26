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

## Migrations

Set:

```bash
export BIOSTAR_DATABASE_URL="postgresql+psycopg://user:password@localhost:5432/biostar"
```

Then run:

```bash
alembic upgrade head
```

The first migration creates the schema only. Scientific reference-data seeding
will be added separately so schema changes and data imports remain distinct.
