from __future__ import annotations

from collections.abc import Iterator

import pytest
from fastapi.testclient import TestClient

from BioStar.engine import AminoAcidData, BiochemistryData
from BioStarAPI.app import app
from BioStarAPI.controllers.analysis import get_sequence_service
from BioStarAPI.services.sequence_service import SequenceService


def _test_data() -> BiochemistryData:
    amino_acids: dict[str, AminoAcidData] = {
        "M": AminoAcidData("M", 149.21, 1.9, 1.2, 0.5, 2.13, 9.28, None),
        "I": AminoAcidData("I", 131.17, 4.5, 1.6, 0.5, 2.36, 9.68, None),
    }
    return BiochemistryData(
        amino_acids=amino_acids,
        aromatic=frozenset(),
        nonpolar=frozenset({"M", "I"}),
        polar=frozenset(),
        positive=frozenset(),
        negative=frozenset(),
        dna_codons={"ATG": "M", "ATT": "I", "ATC": "I", "ATA": "I", "TAA": "*"},
        rna_codons={"AUG": "M", "AUU": "I", "AUC": "I", "AUA": "I", "UAA": "*"},
        dna_stop_codons=frozenset({"TAA"}),
        rna_stop_codons=frozenset({"UAA"}),
        dna_start_codon="ATG",
        rna_start_codon="AUG",
        water_mass=18.01528,
        n_term_pka=7.7,
        c_term_pka=3.5,
    )


class FakeRepository:
    def get_standard_data(self) -> BiochemistryData:
        return _test_data()


@pytest.fixture
def client(monkeypatch: pytest.MonkeyPatch) -> Iterator[TestClient]:
    app.dependency_overrides[get_sequence_service] = lambda: SequenceService(FakeRepository())
    monkeypatch.setattr("BioStarAPI.app.resolve_api_key", lambda raw_key: None)
    monkeypatch.setattr("BioStarAPI.app.enforce_request_limits", lambda request, api_key_id: None)

    with TestClient(app) as test_client:
        yield test_client

    app.dependency_overrides.clear()


def test_dna_to_rna_returns_api_envelope(client: TestClient) -> None:
    response = client.get("/api/v1/dna-rna", params={"sequence": "ATG"})

    assert response.status_code == 200
    assert response.json()["data"] == {"sequence": "AUG"}
    assert response.json()["message"] is None


def test_invalid_api_key_returns_401(client: TestClient) -> None:
    response = client.get(
        "/api/v1/dna-rna",
        params={"sequence": "ATG"},
        headers={"X-API-Key": "invalid"},
    )

    assert response.status_code == 401
    assert response.json()["data"] is None
    assert response.json()["message"]["code"] == "invalid_api_key"


def test_validation_error_preserves_field_details(client: TestClient) -> None:
    response = client.post("/api/v1/protein", json={"get_molecular_weight": True})

    assert response.status_code == 422
    assert response.json()["data"] is None
    message = response.json()["message"]
    assert message["code"] == "validation_error"
    assert message["details"][0]["loc"] == ["body", "sequence"]


def test_sequence_normalization_rejects_invalid_dna(client: TestClient) -> None:
    response = client.get("/api/v1/dna-rna", params={"sequence": "ATX"})

    assert response.status_code == 422
    assert response.json()["message"]["code"] == "invalid_sequence"


def test_dna_to_rna_normalizes_fasta(client: TestClient) -> None:
    response = client.get(
        "/api/v1/dna-rna",
        params={"sequence": ">sequence\nATG"},
    )

    assert response.status_code == 200
    assert response.json()["data"] == {"sequence": "AUG"}


def test_dna_to_rna_rejects_multiple_fasta_records(client: TestClient) -> None:
    response = client.get(
        "/api/v1/dna-rna",
        params={"sequence": ">one\nATG\n>two\nATG"},
    )

    assert response.status_code == 422
    assert response.json()["message"]["code"] == "fasta_single_sequence_required"
