from __future__ import annotations

from collections.abc import Iterator

import pytest
from fastapi.testclient import TestClient

from BioStarAPI.app import app
from BioStarAPI.controllers.analysis import get_sequence_service
from BioStarAPI.services.sequence_service import SequenceService


@pytest.fixture
def client(monkeypatch: pytest.MonkeyPatch) -> Iterator[TestClient]:
    app.dependency_overrides[get_sequence_service] = lambda: SequenceService()
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


def test_batch_dna_to_rna_processes_all_sequences(client: TestClient) -> None:
    response = client.post(
        "/api/v1/batch/dna-rna",
        json={"sequences": ["ATG", "GGC", "TAA"]},
    )

    assert response.status_code == 200
    assert response.json()["data"] == [
        {"sequence": "AUG"},
        {"sequence": "GGC"},
        {"sequence": "UAA"},
    ]


def test_batch_dna_to_rna_rejects_missing_sequences(client: TestClient) -> None:
    response = client.post("/api/v1/batch/dna-rna", json={"sequence": "ATG"})

    assert response.status_code == 422
    assert response.json()["message"]["code"] == "validation_error"


def test_batch_dna_to_rna_rejects_empty_sequences_list(client: TestClient) -> None:
    response = client.post("/api/v1/batch/dna-rna", json={"sequences": []})

    assert response.status_code == 422
    assert response.json()["message"]["code"] == "validation_error"


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
