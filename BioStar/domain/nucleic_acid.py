from __future__ import annotations

import re

from BioStar.data.biochemistry import CODON_SIZE
from BioStar.engine.biochemistry import BiochemistryData
from BioStar.domain.protein import Protein
from BioStar.io.fasta import FastaParserDNA


class NucleicAcid:
    """Base class for DNA and RNA sequence operations."""

    def __init__(self, sequence: str, data: BiochemistryData) -> None:
        self.sequence: str = sequence.upper()
        self.sequence_size: int = len(self.sequence)
        self.sequence_map: dict[str, int] = self.get_sequence_map()
        self.data: BiochemistryData = data

    def get_sequence_map(self) -> dict[str, int]:
        raise NotImplementedError

    def get_peptide_sequence(self, show_stop_codon: bool = False) -> str:
        peptide_sequence: str = ""
        for index in range(0, self.sequence_size - (CODON_SIZE - 1), CODON_SIZE):
            codon: str = self.sequence[index:index + CODON_SIZE]
            peptide_sequence += self.codon_table[codon]

        if show_stop_codon:
            return peptide_sequence
        if peptide_sequence.endswith("*"):
            return peptide_sequence[:-1]
        return peptide_sequence

    def to_protein(self) -> Protein:
        return Protein(self.get_peptide_sequence(), self.data)


class DNA(NucleicAcid):
    """Represents a DNA sequence and DNA-specific analyses."""

    def __init__(self, sequence: str, data: BiochemistryData) -> None:
        normalized_sequence: str = self._fasta_sequence(sequence)
        super().__init__(normalized_sequence, data)
        self.codon_table: dict[str, str] = data.dna_codons

    def get_sequence_map(self) -> dict[str, int]:
        count: dict[str, int] = {"A": 0, "C": 0, "G": 0, "T": 0, "total": 0}
        for nucleotide in self.sequence:
            if nucleotide in count:
                count[nucleotide] += 1
                count["total"] += 1
        return count

    @staticmethod
    def _fasta_sequence(sequence: str = "") -> str:
        return re.sub(r"[^ACGT]", "", sequence.upper())

    def gc_content(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        if self.sequence_size == 0:
            return 0.0
        gc_count: int = self.sequence_map["C"] + self.sequence_map["G"]
        return round((gc_count / self.sequence_size) * multiply_by, decimal_places)

    def at_skew(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        denominator: int = self.sequence_map["A"] + self.sequence_map["T"]
        if denominator == 0:
            return 0.0
        value: float = (self.sequence_map["A"] - self.sequence_map["T"]) / denominator
        return round(value * multiply_by, decimal_places)

    def gc_skew(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        denominator: int = self.sequence_map["G"] + self.sequence_map["C"]
        if denominator == 0:
            return 0.0
        value: float = (self.sequence_map["G"] - self.sequence_map["C"]) / denominator
        return round(value * multiply_by, decimal_places)

    def template_strand(self, reverse_string: bool = True) -> str:
        complement: dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}
        template: str = "".join(complement[nucleotide] for nucleotide in self.sequence)
        if reverse_string:
            return template[::-1]
        return template

    def rna_sequence(self) -> str:
        return self.sequence.replace("T", "U")

    def to_rna(self) -> "RNA":
        return RNA(self.rna_sequence(), self.data)

    def orf_map(self, length_threshold: int = 0) -> list[dict[str, object]]:
        return OpenReadFrame(self.sequence, self.data, length_threshold).orf_map

    def get_orf_map(self, length_threshold: int = 0) -> list[dict[str, object]]:
        return self.orf_map(length_threshold)


class RNA(NucleicAcid):
    """Represents an RNA sequence and RNA-specific analyses."""

    def __init__(self, sequence: str, data: BiochemistryData) -> None:
        normalized_sequence: str = self._fasta_sequence(sequence)
        super().__init__(normalized_sequence, data)
        self.codon_table: dict[str, str] = data.rna_codons

    def get_sequence_map(self) -> dict[str, int]:
        count: dict[str, int] = {"A": 0, "C": 0, "G": 0, "U": 0, "total": 0}
        for nucleotide in self.sequence:
            if nucleotide in count:
                count[nucleotide] += 1
                count["total"] += 1
        return count

    @staticmethod
    def _fasta_sequence(sequence: str = "") -> str:
        return re.sub(r"[^ACGU]", "", sequence.upper())

    def trim_on_stop_codon(self) -> str:
        for index in range(0, self.sequence_size - (CODON_SIZE - 1), CODON_SIZE):
            codon: str = self.sequence[index:index + CODON_SIZE]
            if codon in self.data.rna_stop_codons:
                return self.sequence[:index + CODON_SIZE]
        return self.sequence

    def dna_sequence(self) -> str:
        return self.sequence.replace("U", "T")

    def to_dna(self) -> DNA:
        return DNA(self.dna_sequence(), self.data)


class OpenReadFrame:
    """Finds open reading frames across the six DNA reading frames."""

    def __init__(
        self,
        sequence: str,
        data: BiochemistryData,
        length_threshold: int = 0,
    ) -> None:
        self.length_threshold: int = length_threshold
        self.data: BiochemistryData = data
        self.fasta_map: list[dict[str, str]] = FastaParserDNA().get_sequence_map(sequence)
        self.orf_map: list[dict[str, object]] = []
        self.update_orf_map()
        self.largest_frame: dict[str, int] = self._get_largest_frame()
        self.largest_fame: dict[str, int] = self.largest_frame

    def _get_largest_frame(self) -> dict[str, int]:
        if not self.orf_map:
            return {"nt_length": 0, "aa_length": 0}
        sequence: str = str(self.orf_map[0]["sequence"])
        return {"nt_length": len(sequence), "aa_length": len(sequence) // CODON_SIZE}

    def extract_data_from_fasta_map(
        self,
        sequence: str,
        label: str,
        n: int,
        frame_reference: str,
        sequence_size: int,
    ) -> None:
        frame_sequence: str = ""
        codon_start: int = 0

        for index in range(sequence_size)[n::CODON_SIZE]:
            codon: str = sequence[index:index + CODON_SIZE]
            if not frame_sequence:
                if codon == self.data.dna_start_codon:
                    codon_start = self.find_start_codon(index, frame_reference, sequence_size)
                    frame_sequence = codon
                continue

            frame_sequence += codon
            if codon in self.data.dna_stop_codons:
                if len(frame_sequence) > self.length_threshold:
                    self.orf_map.append({
                        "codon_start": codon_start,
                        "codon_stop": self.find_stop_codon(index, frame_reference, sequence_size),
                        "count_aa": (len(frame_sequence) // CODON_SIZE) - 1,
                        "count_nt": len(frame_sequence) - CODON_SIZE,
                        "label": label,
                        "sequence": frame_sequence,
                        "reference": frame_reference,
                    })
                frame_sequence = ""

    def find_start_codon(self, position: int, frame_reference: str, seq_size: int) -> int:
        if frame_reference in ["+", "++", "+++"]:
            return position + 1
        return abs(position - seq_size)

    def find_stop_codon(self, position: int, frame_reference: str, seq_size: int) -> int:
        if frame_reference in ["+", "++", "+++"]:
            return position + 2
        return abs((position + 2) - seq_size)

    def update_orf_map(self) -> list[dict[str, object]]:
        frame_references: list[str] = ["+", "++", "+++", "-", "--", "---"]

        for item in self.fasta_map:
            dna: DNA = DNA(item["sequence"], self.data)
            for frame_index in range(CODON_SIZE):
                self.extract_data_from_fasta_map(
                    dna.sequence, item["label"], frame_index,
                    frame_references[frame_index], dna.sequence_size,
                )
            template: str = dna.template_strand()
            for frame_index in range(CODON_SIZE):
                self.extract_data_from_fasta_map(
                    template, item["label"], frame_index,
                    frame_references[frame_index + CODON_SIZE], dna.sequence_size,
                )

        self.orf_map = sorted(
            self.orf_map,
            key=lambda item: len(str(item["sequence"])),
            reverse=True,
        )
        return self.orf_map
