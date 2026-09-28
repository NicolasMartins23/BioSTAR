from __future__ import annotations

from BioStar.Protein.protein import Protein
from BioStarAPI.controllers.schemas import ProteinAnalysisRequest
from BioStarAPI.database.repositories.biochemistry import BiochemistryRepository
from BioStarAPI.services.sequence_service import SequenceService


class ProteinService:
    def __init__(self, repository: BiochemistryRepository) -> None:
        self.repository = repository
        self.sequence_service = SequenceService(repository)

    def analyze(self, request: ProteinAnalysisRequest) -> dict[str, object]:
        data = self.repository.get_standard_data()
        sequence = self.sequence_service.normalize_protein(
            request.sequence,
            10_000,
        )
        protein = Protein(sequence, data)

        requested = {
            "aminoacids_count": request.get_aminoacids_count,
            "isoelectric_point": request.get_isoelectric_point,
            "charge_at_pH": request.get_charge_at_pH is not None,
            "aromaticity": request.get_aromaticity,
            "secondary_structure_propensity": request.get_secondary_structure_propensity,
            "molecular_weight": request.get_molecular_weight,
            "hydrophobic_index": request.get_hydrophobic_index,
            "composition_ratio": request.get_composition_ratio,
            "extinction_coefficient": request.get_extinction_coefficient,
        }
        if request.get_full_test_results:
            requested = {name: True for name in requested}

        results: dict[str, object] = {"sequence": protein.sequence, "length": protein.sequence_size}
        if requested["aminoacids_count"]:
            results["aminoacids_count"] = protein.count["by_aminoacid"]
        if requested["isoelectric_point"]:
            results["isoelectric_point"] = protein.isoelectric_point()
        if requested["charge_at_pH"]:
            pH = request.get_charge_at_pH if request.get_charge_at_pH is not None else 7.0
            results["charge_at_pH"] = {"pH": pH, "charge": protein.charge_at_pH(pH)}
        if requested["aromaticity"]:
            results["aromaticity"] = protein.aromacity()
        if requested["secondary_structure_propensity"]:
            results["secondary_structure_propensity"] = protein.secondary_structure_propensity()
        if requested["molecular_weight"]:
            results["molecular_weight"] = protein.molecular_weight()
        if requested["hydrophobic_index"]:
            results["hydrophobic_index"] = protein.hydrophobic_index()
        if requested["composition_ratio"]:
            results["composition_ratio"] = protein.composition_ratio()
        if requested["extinction_coefficient"]:
            results["extinction_coefficient"] = protein.extinction_coefficient()
        return results
