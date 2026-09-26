from re import sub

from BioStar.engine.biochemistry import BiochemistryData
from BioStar.engine.default_data import get_default_biochemistry

WATER_MASS = 18.015


class Protein:
    """Represents a protein sequence and provides sequence analyses."""

    def __init__(self, sequence: str = "", data: BiochemistryData | None = None) -> None:
        self.data: BiochemistryData = data or get_default_biochemistry()
        self.sequence: str = sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "", sequence.upper())
        self.count: dict[str, int | dict[str, int]] = self._update_count()
        self.sequence_size: int = int(self.count["total"])

    def _update_count(self) -> dict[str, int | dict[str, int]]:
        count: dict[str, int] = {aa: 0 for aa in self.data.amino_acids if aa != "*"}
        aromatic = nonpolar = polar = polar_negative = polar_neutral = polar_positive = 0

        for aa in self.sequence:
            count[aa] += 1
            if aa in self.data.aromatic:
                aromatic += 1
            elif aa in self.data.nonpolar:
                nonpolar += 1
            elif aa in self.data.polar:
                polar += 1
                if aa in self.data.negative:
                    polar_negative += 1
                elif aa in self.data.positive:
                    polar_positive += 1
                else:
                    polar_neutral += 1

        return {
            "by_aminoacid": count,
            "aromatic": aromatic,
            "nonpolar": nonpolar,
            "polar": polar,
            "polar_negative": polar_negative,
            "polar_neutral": polar_neutral,
            "polar_positive": polar_positive,
            "total": len(self.sequence),
        }

    def _aminoacid_counts(self) -> dict[str, int]:
        return self.count["by_aminoacid"]  # type: ignore[return-value]

    def aromacity(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        if self.sequence_size == 0:
            return 0.0
        return round((int(self.count["aromatic"]) / self.sequence_size) * multiply_by, decimal_places)

    def charge_at_pH(self, pH: float = 7.0) -> float:
        normalized_pH = round(pH, 3)
        positive = 0.0
        negative = 0.0
        counts = self._aminoacid_counts()

        for aa in self.data.positive:
            pKr = self.data.amino_acids[aa].pkr
            if pKr is not None:
                positive += counts[aa] / (1.0 + 10 ** (normalized_pH - pKr))

        for aa in self.data.negative:
            pKr = self.data.amino_acids[aa].pkr
            if pKr is not None:
                negative += counts[aa] / (1.0 + 10 ** (pKr - normalized_pH))

        positive += 1.0 / (1.0 + 10 ** (normalized_pH - 9.69))
        negative += 1.0 / (1.0 + 10 ** (2.34 - normalized_pH))
        return round(positive - negative, 2)

    def composition_ratio(self, multiply_by: float = 1.0, decimal_places: int = 4) -> dict[str, float]:
        counts = self._aminoacid_counts()
        if self.sequence_size == 0:
            return {aa: 0.0 for aa in counts}
        return {
            aa: round((count / self.sequence_size) * multiply_by, decimal_places)
            for aa, count in counts.items()
        }

    def extinction_coefficient(self) -> dict[str, float]:
        cys = self.sequence.count("C")
        coefficient = (self.sequence.count("W") * 5500) + (self.sequence.count("Y") * 1490)
        return {
            "cys_cystines": round(coefficient + ((cys // 2) * 125), 2),
            "cys_reduced": round(coefficient, 2),
        }

    def hydrophobic_index(self) -> float:
        if self.sequence_size == 0:
            return 0.0
        total = sum(self.data.amino_acids[aa].hydrophobicity for aa in self.sequence)
        return round(total / self.sequence_size, 2)

    def isoelectric_point(self) -> float:
        low = 0.0
        high = 14.0

        while high - low > 0.01:
            middle = (low + high) / 2.0
            charge = self.charge_at_pH(middle)
            if abs(charge) < 0.01:
                return round(middle, 2)
            if charge > 0:
                low = middle
            else:
                high = middle

        return round((low + high) / 2.0, 2)

    def molecular_weight(self) -> float:
        if self.sequence_size == 0:
            return 0.0
        weight = sum(self.data.amino_acids[aa].molecular_weight for aa in self.sequence)
        weight -= (self.sequence_size - 1) * WATER_MASS
        return round(weight, 2)

    def secondary_structure_propensity(self) -> dict[str, float]:
        if self.sequence_size == 0:
            return {"alpha_helix": 0.0, "beta_sheet": 0.0, "coil": 0.0}

        alpha = sum(self.data.amino_acids[aa].alpha_helix for aa in self.sequence) / self.sequence_size
        beta = sum(self.data.amino_acids[aa].beta_sheet for aa in self.sequence) / self.sequence_size
        return {
            "alpha_helix": round(alpha * 100, 2),
            "beta_sheet": round(beta * 100, 2),
            "coil": round((1 - alpha - beta) * 100, 2),
        }
