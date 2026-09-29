from __future__ import annotations

import re

from BioStarAPI.resources.exceptions import BioStarAPIError
from BioStarAPI.resources.messages import MessageCode

from BioStar.NucleicAcids.nucleic_acids import DNA, RNA


DNA_ALPHABET = frozenset("ACGT")
RNA_ALPHABET = frozenset("ACGU")
PROTEIN_ALPHABET = frozenset("ACDEFGHIKLMNPQRSTVWY")


class SequenceService:
    def normalize(self, sequence: str, alphabet: frozenset[str], kind: str, max_length: int) -> str:
        value = sequence.strip().upper()
        if value.startswith(">"):
            records = self._parse_fasta(value)
            if len(records) != 1:
                raise BioStarAPIError(422, MessageCode.FASTA_SINGLE_SEQUENCE_REQUIRED)
            value = records[0]
        else:
            value = re.sub(r"\s+", "", value)

        if not value:
            raise BioStarAPIError(422, MessageCode.EMPTY_SEQUENCE, kind=kind)
        if len(value) > max_length:
            raise BioStarAPIError(422, MessageCode.SEQUENCE_TOO_LONG, kind=kind, max_length=max_length)
        invalid = sorted(set(value) - alphabet)
        if invalid:
            raise BioStarAPIError(422, MessageCode.INVALID_SEQUENCE, kind=kind, characters=", ".join(invalid))
        return value

    def normalize_dna(self, sequence: str, max_length: int) -> str:
        return self.normalize(sequence, DNA_ALPHABET, "DNA", max_length)

    def normalize_rna(self, sequence: str, max_length: int) -> str:
        return self.normalize(sequence, RNA_ALPHABET, "RNA", max_length)

    def normalize_protein(self, sequence: str, max_length: int) -> str:
        return self.normalize(sequence, PROTEIN_ALPHABET, "protein", max_length)

    def dna_to_rna(self, sequence: str, max_length: int) -> dict[str, str]:
        normalized = self.normalize_dna(sequence, max_length)
        return {"sequence": DNA(normalized).rna_sequence()}

    def dna_to_protein(self, sequence: str, max_length: int) -> dict[str, str]:
        normalized = self.normalize_dna(sequence, max_length)
        return {"sequence": DNA(normalized).to_protein().sequence}

    def rna_to_protein(self, sequence: str, max_length: int) -> dict[str, str]:
        normalized = self.normalize_rna(sequence, max_length)
        return {"sequence": RNA(normalized).to_protein().sequence}

    def rna_to_dna(self, sequence: str, max_length: int) -> dict[str, str]:
        normalized = self.normalize_rna(sequence, max_length)
        return {"sequence": RNA(normalized).dna_sequence()}

    @staticmethod
    def _parse_fasta(value: str) -> list[str]:
        records: list[str] = []
        current: list[str] | None = None
        for line in value.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    records.append("".join(current))
                current = []
                continue
            if current is None:
                raise BioStarAPIError(422, MessageCode.INVALID_FASTA)
            current.append(line)
        if current is not None:
            records.append("".join(current))
        return records
