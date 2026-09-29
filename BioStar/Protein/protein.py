from __future__ import annotations

from re import sub

from BioStar.engine import get_reference_data


class Protein:
    """Represents a protein sequence and provides sequence analyses."""

    def __init__(self, sequence: str = "") -> None:
        self.reference_data = get_reference_data()
        self.sequence: str = sub(
            r"[^ACDEFGHIKLMNPQRSTVWY]",
            "",
            sequence.upper(),
        )
        self.count: dict[str, int | dict[str, int]] = self._update_count()
        self.sequence_size: int = int(self.count["total"])

    def _update_count(self) -> dict[str, int | dict[str, int]]:
        count: dict[str, int] = {aa: 0 for aa in self.reference_data.amino_acid_symbols()}
        aromatic: int = 0
        nonpolar: int = 0
        polar: int = 0
        polar_negative: int = 0
        polar_neutral: int = 0
        polar_positive: int = 0

        for aa in self.sequence:
            count[aa] += 1
            if aa in self.reference_data.amino_acid_class("aromatic"):
                aromatic += 1
            elif aa in self.reference_data.amino_acid_class("nonpolar"):
                nonpolar += 1
            elif aa in self.reference_data.amino_acid_class("polar"):
                polar += 1
                if aa in self.reference_data.amino_acid_class("negative"):
                    polar_negative += 1
                elif aa not in self.reference_data.amino_acid_class("positive"):
                    polar_neutral += 1
                else:
                    polar_positive += 1

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
        normalized_pH: float = round(pH, 3)
        positive: float = 0.0
        negative: float = 0.0
        counts: dict[str, int] = self._aminoacid_counts()

        for aa in self.reference_data.amino_acid_class("positive"):
            pKa = self.reference_data.amino_acid(aa).pkr
            if pKa is not None:
                positive += counts[aa] / (1.0 + 10 ** (normalized_pH - pKa))

        for aa in self.reference_data.amino_acid_class("negative"):
            pKa = self.reference_data.amino_acid(aa).pkr
            if pKa is not None:
                negative += counts[aa] / (1.0 + 10 ** (pKa - normalized_pH))

        positive += 1.0 / (1.0 + 10 ** (normalized_pH - self.reference_data.constant("n_term_pka")))
        negative += 1.0 / (1.0 + 10 ** (self.reference_data.constant("c_term_pka") - normalized_pH))
        return round(positive - negative, 2)

    def composition_ratio(self, multiply_by: float = 1.0, decimal_places: int = 4) -> dict[str, float]:
        counts: dict[str, int] = self._aminoacid_counts()
        if self.sequence_size == 0:
            return {aa: 0.0 for aa in counts}
        return {
            aa: round((count / self.sequence_size) * multiply_by, decimal_places)
            for aa, count in counts.items()
        }

    def extinction_coefficient(self) -> dict[str, float]:
        cys: int = self.sequence.count("C")
        coefficient: float = (self.sequence.count("W") * 5500) + (self.sequence.count("Y") * 1490)
        return {
            "cys_cystines": round(coefficient + ((cys // 2) * 125), 2),
            "cys_reduced": round(coefficient, 2),
        }

    def hydrophobic_index(self) -> float:
        if self.sequence_size == 0:
            return 0.0
        total: float = sum(
            self.reference_data.amino_acid(aa).hydrophobicity
            for aa in self.sequence
        )
        return round(total / self.sequence_size, 2)

    def isoelectric_point(self) -> float:
        low: float = 0.0
        high: float = 14.0

        while high - low > 0.01:
            middle: float = (low + high) / 2.0
            charge: float = self.charge_at_pH(middle)
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
        weight: float = sum(
            self.reference_data.amino_acid(aa).molecular_weight
            for aa in self.sequence
        )
        weight -= (self.sequence_size - 1) * self.reference_data.constant("water_mass")
        return round(weight, 2)

    def secondary_structure_propensity(self) -> dict[str, float]:
        if self.sequence_size == 0:
            return {"alpha_helix": 0.0, "beta_sheet": 0.0, "coil": 0.0}

        alpha: float = sum(
            self.reference_data.amino_acid(aa).alpha_helix
            for aa in self.sequence
        ) / self.sequence_size
        beta: float = sum(
            self.reference_data.amino_acid(aa).beta_sheet
            for aa in self.sequence
        ) / self.sequence_size
        return {
            "alpha_helix": round(alpha * 100, 2),
            "beta_sheet": round(beta * 100, 2),
            "coil": round((1 - alpha - beta) * 100, 2),
        }
