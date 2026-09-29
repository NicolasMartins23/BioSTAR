from __future__ import annotations

import os
from collections.abc import Generator

from sqlalchemy import create_engine
from sqlalchemy.orm import Session


DATABASE_URL: str = os.getenv(
    "BIOSTAR_DATABASE_URL",
    "postgresql+psycopg://biostar:biostar@db:5432/biostar",
)

engine = create_engine(DATABASE_URL, pool_pre_ping=True)


def get_session() -> Generator[Session, None, None]:
    with Session(engine) as session:
        yield session
