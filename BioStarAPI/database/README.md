# BioSTAR API Database

The PostgreSQL database belongs to the API application. The BioSTAR engine
does not depend on SQLAlchemy, PostgreSQL, or FastAPI.

## PostgreSQL data

The API database stores application state such as:

- users
- refresh tokens
- API keys
- API request usage

Scientific reference data is owned by the BioSTAR engine and bundled as a
SQLite database under `BioStar/data/`.

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

No PostgreSQL seed is required for biochemical reference data.
