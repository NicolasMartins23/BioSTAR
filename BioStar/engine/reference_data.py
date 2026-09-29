from __future__ import annotations

from functools import lru_cache

from BioStar.engine.biochemistry import AminoAcidData
from BioStar.engine.default_data import get_default_biochemistry


class ReferenceData:
    """Provides read-only access to bundled biochemical reference data."""

    def __init__(self) -> None:
        self._data = get_default_biochemistry()

    def amino_acid(self, symbol: str) -> AminoAcidData:
        try:
            return self._data.amino_acids[symbol]
        except KeyError:
            raise KeyError(f"Unknown amino acid: {symbol}.") from None

    def amino_acid_symbols(self) -> tuple[str, ...]:
        return tuple(sorted(self._data.amino_acids))

    def amino_acid_class(self, class_name: str) -> frozenset[str]:
        classes = {
            "aromatic": self._data.aromatic,
            "nonpolar": self._data.nonpolar,
            "polar": self._data.polar,
            "positive": self._data.positive,
            "negative": self._data.negative,
        }
        try:
            return classes[class_name]
        except KeyError:
            raise KeyError(f"Unknown amino acid class: {class_name}.") from None

    def codon_table(self, molecule_type: str) -> dict[str, str]:
        tables = {
            "DNA": self._data.dna_codons,
            "RNA": self._data.rna_codons,
        }
        try:
            return dict(tables[molecule_type])
        except KeyError:
            raise KeyError(f"Unknown molecule type: {molecule_type}.") from None

    def stop_codons(self, molecule_type: str) -> frozenset[str]:
        tables = {
            "DNA": self._data.dna_stop_codons,
            "RNA": self._data.rna_stop_codons,
        }
        try:
            return tables[molecule_type]
        except KeyError:
            raise KeyError(f"Unknown molecule type: {molecule_type}.") from None

    def start_codon(self, molecule_type: str) -> str:
        codons = {
            "DNA": self._data.dna_start_codon,
            "RNA": self._data.rna_start_codon,
        }
        try:
            return codons[molecule_type]
        except KeyError:
            raise KeyError(f"Unknown molecule type: {molecule_type}.") from None

    def constant(self, key: str) -> float:
        constants = {
            "water_mass": self._data.water_mass,
            "n_term_pka": self._data.n_term_pka,
            "c_term_pka": self._data.c_term_pka,
        }
        try:
            return float(constants[key])
        except KeyError:
            raise KeyError(f"Unknown biochemical reference constant: {key}.") from None

    def close(self) -> None:
        """Retained as a no-op for backwards compatibility."""
        return None

    def __enter__(self) -> "ReferenceData":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


@lru_cache(maxsize=1)
def get_reference_data() -> ReferenceData:
    return ReferenceData()
