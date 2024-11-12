import re
from .DataBiochemistry import (
    CODON_SIZE,
    TABLE_DNA_CODON_TO_AMINOACID, STOP_CODON_DNA,
    TABLE_RNA_CODON_TO_AMINOACID, STOP_CODON_RNA,
    )

from .protein import Protein
from .tools.orf import OpenReadFrame
from .data_handler import FastaParserDNA

class NucleicAcid():
    """
    Base class for Nucleic Acid functionality.
    The DNA and RNA classes use this as a parent class.
    This class has functions that handle the logic of biological processes DNA and RNA share.
    """
    def __init__(
            self,
            sequence: str,
            codon_table: dict,
        ) -> None:
        self.sequence = sequence
        self.codon_table = codon_table
        self.count = self.get_count()
        self.sequence_size = self.count["total"]


    def get_count(self) -> dict:
        """
        This function will be overritten.
        It must return a dictionary with the count details
        A count, C count, G count, T/U count,
        """
        pass


    def get_peptide_sequence(
            self,
            show_stop_codon: bool = False
        ) -> str:
        """
        Returns the peptide sequence string of a NucleicAcid instance.
        The codon_table is a variable that is defined on the children class callback of the function.
        If you wish to return a protein object instead, use the method to_protein()
        """
        peptide_sequence: str = ""

        for i in range(0, self.sequence_size, CODON_SIZE):
            codon: str = self.sequence[i:i + CODON_SIZE]
            peptide_sequence += self.codon_table[codon]

        if show_stop_codon:
            return peptide_sequence
        else:
            return peptide_sequence[:-1]
        
    def to_protein(self) -> Protein:
        peptide_sequence = self.get_peptide_sequence()
        return Protein(peptide_sequence)


class DNA(NucleicAcid):
    """
    This class handles DNA specialized functions and inherits NucleicAcid functions.
    It must be initialized passing a valid nucleotide sequence as a parameter.
    """
    def __init__(
            self,
            sequence: str,
        ) -> None:
        super().__init__(sequence, TABLE_DNA_CODON_TO_AMINOACID)

    def get_count(self) -> dict:
        # Overrides parent function
        count: dict = {
            "A": 0,
            "C": 0,
            "G": 0,
            "T": 0,
            "total": 0
        }
        for nucleotide in self.sequence:
            count["total"] += 1
            if nucleotide == "A":
                count["A"] += 1
            elif nucleotide == "C":
                count["C"] += 1
            elif nucleotide == "G":
                count["G"] += 1
            elif nucleotide == "T":
                count["T"] += 1
        return count

    def _fasta_sequence(self, dna_seq: str = "") -> str:
        fasta_sequence = re.sub(r'[^ACGT]', '', dna_seq.upper())
        return fasta_sequence


    def get_gc_content(
            self,
            multiply_result_by_100: bool = True,
            decimal_digits: int = 1
            ) -> float:
        """
        When multiply_result_by_100 is true the result returns as 40, 60, etc.
        Else, it will return as 0.4 0.6
        """
        gc_count: int = self.count["C"] + self.count["G"]
        gc_content: float = (gc_count / self.sequence_size)

        if multiply_result_by_100:            
            return round((gc_content * 100), decimal_digits)
        else:
            return round((gc_content), decimal_digits)

    def get_at_skew(
            self,
            multiply_result_by_100: bool = True,
            decimal_digits: int = 1
            ) -> float:
        """
        The GC skew is positive and negative in the leading strand and in the lagging strand respectively; therefore, it is expected to see a switch in GC skew sign just at the point of DNA replication origin and terminus.
        When multiply_result_by_100 is true the result returns as 40, 60, etc.
        Else, it will return as 0.4 0.6
        """
        at_skew: float = (self.count["A"] - self.count["T"]) / (self.count["A"] + self.count["T"])

        if multiply_result_by_100:            
            return round((at_skew * 100), decimal_digits)
        else:
            return round((at_skew), decimal_digits)

    def get_gc_skew(
            self,
            multiply_result_by_100: bool = True,
            decimal_digits: int = 1
            ) -> float:
        """
        The GC skew is positive and negative in the leading strand and in the lagging strand respectively; therefore, it is expected to see a switch in GC skew sign just at the point of DNA replication origin and terminus.
        When multiply_result_by_100 is true the result returns as 40, 60, etc.
        Else, it will return as 0.4 0.6
        """
        gc_skew: float = (self.count["G"] - self.count["C"]) / (self.count["G"] + self.count["C"])

        if multiply_result_by_100:            
            return round((gc_skew * 100), decimal_digits)
        else:
            return round((gc_skew), decimal_digits)

    def get_template_strand(self,
                            reverse_string: bool = True) -> str:
        template_strand: str = ""
        for nucleotide in self.sequence:
            if nucleotide == "A":
                template_strand += "T"
            elif nucleotide == "T":
                template_strand += "A"
            elif nucleotide == "C":
                template_strand += "G"
            elif nucleotide == "G":
                template_strand += "C"

        if reverse_string:
            return template_strand[::-1]
        else:
            return template_strand

    def get_rna_string(self) -> str:
        """
        Returns a DNA sequence of a given RNA instance.\n
        """
        rna_sequence: str = self.sequence.replace("T", "U")
        return rna_sequence

    def to_rna(self):
        rna_sequence: str = self.get_rna_string()
        return RNA(rna_sequence)

    def get_orf_map(self, length_threshold: int = 0) -> dict:
        orf = OpenReadFrame(sequence=self.sequence, length_threshold=length_threshold)
        return orf.orf_map


class RNA(NucleicAcid):
    def __init__(
            self,
            sequence: str,
        ) -> None:
        super().__init__(sequence, TABLE_RNA_CODON_TO_AMINOACID)

    def get_count(self) -> dict:
        # Overrides parent function
        count: dict = {
            "A": 0,
            "C": 0,
            "G": 0,
            "U": 0,
            "total": 0
        }
        for nucleotide in self.sequence:
            count["total"] += 1
            if nucleotide == "A":
                count["A"] += 1
            elif nucleotide == "C":
                count["C"] += 1
            elif nucleotide == "G":
                count["G"] += 1
            elif nucleotide == "U":
                count["U"] += 1
        return count

    def _fasta_sequence(self, sequence: str = "") -> str:
        fasta_sequence = re.sub(r'[^ACGU]', '', sequence.upper())
        return fasta_sequence

    def trim_on_stop_codon(self):
        trimmed_sequence: str = ""
        for i in range(0, self.sequence_size, CODON_SIZE):
                codon: str = self.sequence[i:i + CODON_SIZE]
                trimmed_sequence += codon
                if codon in STOP_CODON_RNA:
                    return trimmed_sequence

    def get_dna_string(self) -> str:
        """
        Returns a DNA sequence of a given RNA instance.\n
        """
        dna_sequence: str = self.sequence.replace("U", "T")
        return dna_sequence

    def to_dna(self):
        dna_sequence: str = self.get_dna_string()
        return DNA(dna_sequence)
