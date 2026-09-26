from __future__ import annotations

from fastapi import Depends
from sqlalchemy.orm import Session

from BioStarAPI.database.connection import get_session
from BioStarAPI.database.repositories.biochemistry import BiochemistryRepository


def get_biochemistry_repository(
    session: Session = Depends(get_session),
) -> BiochemistryRepository:
    return BiochemistryRepository(session)
