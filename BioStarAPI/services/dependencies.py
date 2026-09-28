from __future__ import annotations

from BioStarAPI.services.sequence_service import SequenceService


def get_sequence_service() -> SequenceService:
    return SequenceService()
