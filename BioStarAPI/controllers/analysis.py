from __future__ import annotations

from fastapi import Depends

from BioStarAPI.controllers.schemas import BatchSequenceRequest, MutationCompareRequest, ProteinAnalysisRequest
from BioStarAPI.database.repositories.biochemistry import BiochemistryRepository
from BioStarAPI.services.dependencies import get_biochemistry_repository
from BioStarAPI.services.mutation_service import MutationService
from BioStarAPI.services.protein_service import ProteinService
from BioStarAPI.services.sequence_service import SequenceService

GET_SEQUENCE_MAX_LENGTH = 1_000
POST_SEQUENCE_MAX_LENGTH = 10_000


def get_sequence_service(
    repository: BiochemistryRepository = Depends(get_biochemistry_repository),
) -> SequenceService:
    return SequenceService(repository)


def dna_to_rna(
    sequence: str,
    max_length: int = POST_SEQUENCE_MAX_LENGTH,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.dna_to_rna(sequence, max_length)


def dna_to_protein(
    sequence: str,
    max_length: int = POST_SEQUENCE_MAX_LENGTH,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.dna_to_protein(sequence, max_length)


def rna_to_protein(
    sequence: str,
    max_length: int = POST_SEQUENCE_MAX_LENGTH,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.rna_to_protein(sequence, max_length)


def rna_to_dna(
    sequence: str,
    max_length: int = POST_SEQUENCE_MAX_LENGTH,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.rna_to_dna(sequence, max_length)


def protein_analysis(
    request: ProteinAnalysisRequest,
    repository: BiochemistryRepository = Depends(get_biochemistry_repository),
) -> dict[str, object]:
    return ProteinService(repository).analyze(request)


def mutation_compare(
    reference: str,
    sequence: str,
    repository: BiochemistryRepository = Depends(get_biochemistry_repository),
) -> dict[str, object]:
    return MutationService(repository).compare(reference, sequence)


def batch_dna_rna(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.dna_to_rna(request.sequence, POST_SEQUENCE_MAX_LENGTH)


def batch_dna_protein(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.dna_to_protein(request.sequence, POST_SEQUENCE_MAX_LENGTH)


def batch_rna_protein(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.rna_to_protein(request.sequence, POST_SEQUENCE_MAX_LENGTH)


def batch_rna_dna(
    request: BatchSequenceRequest,
    service: SequenceService = Depends(get_sequence_service),
) -> dict[str, str]:
    return service.rna_to_dna(request.sequence, POST_SEQUENCE_MAX_LENGTH)
