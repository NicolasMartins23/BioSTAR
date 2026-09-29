from BioStarAPI.controllers.schemas import ProteinAnalysisRequest
from BioStarAPI.services.mutation_service import MutationService
from BioStarAPI.services.protein_service import ProteinService
from BioStarAPI.services.sequence_service import SequenceService


def test_sequence_service_converts_dna_to_rna() -> None:
    service = SequenceService()

    assert service.dna_to_rna("ATG", 1_000) == {"sequence": "AUG"}


def test_sequence_service_converts_dna_to_protein() -> None:
    service = SequenceService()

    assert service.dna_to_protein("ATGGGC", 1_000) == {"sequence": "MG"}


def test_sequence_service_converts_rna_to_protein() -> None:
    service = SequenceService()

    assert service.rna_to_protein("AUGGGC", 1_000) == {"sequence": "MG"}


def test_sequence_service_converts_rna_to_dna() -> None:
    service = SequenceService()

    assert service.rna_to_dna("AUG", 1_000) == {"sequence": "ATG"}


def test_protein_service_analyzes_requested_result() -> None:
    service = ProteinService()
    request = ProteinAnalysisRequest(
        sequence="MG",
        get_molecular_weight=True,
    )

    result = service.analyze(request)

    assert result["sequence"] == "MG"
    assert result["length"] == 2
    assert result["molecular_weight"] == 206.27


def test_mutation_service_returns_only_mutations() -> None:
    service = MutationService()

    result = service.compare("ATG", "AAG")

    assert result["reference"] == "ATG"
    assert result["sequence"] == "AAG"
    assert len(result["mutations"]) == 1
    assert result["mutations"][0]["mutation_type"] == "missense"
