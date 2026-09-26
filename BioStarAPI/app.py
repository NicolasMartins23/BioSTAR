from __future__ import annotations

from fastapi import FastAPI
from sqlalchemy import text

from BioStarAPI.database.connection import engine

app = FastAPI(title="BioSTAR API", version="0.2.0")


@app.get("/health")
def health() -> dict[str, str]:
    with engine.connect() as connection:
        connection.execute(text("SELECT 1"))
    return {"status": "ok"}
