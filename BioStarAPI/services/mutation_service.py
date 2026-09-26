from __future__ import annotations

from fastapi import HTTPException

from BioStar.analysis.sequence_comparison import CompareNucleotideSequence
from BioStarAPI.database.repositories.biochemistry import BiochemistryRepository
from BioStarAPI.services.sequence_service import SequenceService


class MutationService:
    def __init__(self, repository: BiochemistryRepository) -> None:
        self.sequence_service = SequenceService(repository)
        self.repository = repository

    def compare(self, reference: str, sequence: str) -> dict[str, object]:
        reference_dna = self.sequence_service.normalize_dna(reference, 10_000)
        sequence_dna = self.sequence_service.normalize_dna(sequence, 10_000)

        if len(reference_dna) != len(sequence_dna):
            raise HTTPException(
                status_code=422,
                detail="Reference and sequence must have the same length",
            )
        if len(reference_dna) % 3 != 0:
            raise HTTPException(
                status_code=422,
                detail="Reference and sequence lengths must be multiples of 3",
            )

        data = self.repository.get_standard_data()
        mutations = CompareNucleotideSequence(
            reference_dna,
            sequence_dna,
            data,
        ).compare(show_only_mutations=True)

        return {
            "reference": reference_dna,
            "sequence": sequence_dna,
            "mutations": mutations,
        }
