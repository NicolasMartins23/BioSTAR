from __future__ import annotations

from pydantic import BaseModel, Field, StrictBool


class ProteinAnalysisRequest(BaseModel):
    sequence: str = Field(min_length=1)
    get_full_test_results: StrictBool = False
    get_aminoacids_count: StrictBool = False
    get_isoelectric_point: StrictBool = False
    get_charge_at_pH: float | None = None
    get_aromaticity: StrictBool = False
    get_secondary_structure_propensity: StrictBool = False
    get_molecular_weight: StrictBool = False
    get_hydrophobic_index: StrictBool = False
    get_composition_ratio: StrictBool = False
    get_extinction_coefficient: StrictBool = False


class MutationCompareRequest(BaseModel):
    reference: str = Field(min_length=1)
    sequence: str = Field(min_length=1)


class BatchSequenceRequest(BaseModel):
    sequence: str = Field(min_length=1)
