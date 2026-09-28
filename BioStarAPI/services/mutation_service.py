from __future__ import annotations

from BioStar.analysis.sequence_comparison import CompareNucleotideSequence
from BioStar.engine import get_default_biochemistry
from BioStarAPI.resources.exceptions import BioStarAPIError
from BioStarAPI.resources.messages import MessageCode
from BioStarAPI.services.sequence_service import SequenceService


class MutationService:
    def __init__(self) -> None:
        self.sequence_service = SequenceService()

    def compare(self, reference: str, sequence: str) -> dict[str, object]:
        reference_dna = self.sequence_service.normalize_dna(reference, 10_000)
        sequence_dna = self.sequence_service.normalize_dna(sequence, 10_000)

        if len(reference_dna) != len(sequence_dna):
            raise BioStarAPIError(
                422,
                MessageCode.MUTATION_SEQUENCES_SAME_LENGTH,
            )
        if len(reference_dna) % 3 != 0:
            raise BioStarAPIError(
                422,
                MessageCode.MUTATION_SEQUENCES_MULTIPLE_OF_THREE,
            )

        data = get_default_biochemistry()
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
