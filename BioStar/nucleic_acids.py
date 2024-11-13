import re
from .DataBiochemistry import (CODON_SIZE,
    TABLE_DNA_CODON_TO_AMINOACID, STOP_CODON_DNA,
    TABLE_RNA_CODON_TO_AMINOACID, STOP_CODON_RNA,
    )

from .protein import Protein

from .data_handler import FastaParserDNA

class NucleicAcid():
    """
    Base class for Nucleic Acid functionality.
    The DNA and RNA classes use this as a parent class.
    This class has functions that handle the logic of biological processes DNA and RNA share.
    """
    def __init__(self, sequence: str, codon_table: dict) -> None:
        self.sequence = sequence
        self.sequence_size = len(self.sequence)
        self.sequence_map = self.get_sequence_map()
        self.codon_table = codon_table

    def get_peptide_sequence(self, show_stop_codon: bool = False) -> str:
        """
        Returns a string which represents the peptide sequence of a given NucleicAcid when translated intto a protein.
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
        """
        Returns a protein object of a given NucleicAcid when translated into a protein.
        If you wish to return a protein object instead, use the method to_protein()
        """
        peptide_sequence = self.get_peptide_sequence()
        return Protein(peptide_sequence)


class DNA(NucleicAcid):
    """
    This class handles DNA specific functions and inherits NucleicAcid functions.
    It must be initialized passing a valid nucleotide sequence as a parameter.
    """
    def __init__(self, sequence: str) -> None:
        super().__init__(sequence, TABLE_DNA_CODON_TO_AMINOACID)

    def get_sequence_map(self) -> dict:
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


    def gc_content(self, multiply_by: float = 1.0) -> float:
        """
        Returns the ratio of Guanine and Cytosine content of the sequence.
        """
        gc_count: int = self.count["C"] + self.count["G"]
        gc_content: float = (gc_count / self.sequence_size) * multiply_by

        return round(gc_content, 2)

    def at_skew(self, multiply_by: float = 1.0) -> float:
        """
        Returns the skew of Adenine and Thymine content of the sequence.
        A higher value means a prevalence of Adenine.
        A lower value means a prevalence of Cytosine.
        """
        at_skew: float = (self.count["A"] - self.count["T"]) / (self.count["A"] + self.count["T"])

        return round((at_skew * multiply_by), 2)

    def gc_skew(self, multiply_by: float = 1.0) -> float:
        """
        Returns the skew of Guanine and Cytosine content of the sequence.
        A higher value means a prevalence of Guanine.
        A lower value means a prevalence of Cytosine.
        """
        gc_skew: float = (self.count["G"] - self.count["C"]) / (self.count["G"] + self.count["C"])

        if multiply_result_by_100:            
            return round((gc_skew * 100), decimal_digits)
        else:
            return round((gc_skew), decimal_digits)

    def template_strand(self, reverse_string: bool = True) -> str:
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

    def rna_sequence(self) -> str:
        """
        Returns a string which represents the corresponding RNA sequence.
        """
        rna_sequence: str = self.sequence.replace("T", "U")
        return rna_sequence

    def to_rna(self):
        """
        Returns an RNA() object instance.
        """
        rna_sequence: str = self.rna_sequence()
        return RNA(rna_sequence)

    def orf_map(self, length_threshold: int = 0) -> dict:
        orf = OpenReadFrame(sequence=self.sequence, length_threshold=length_threshold)
        return orf.orf_map


class RNA(NucleicAcid):
    def __init__(self, sequence: str) -> None:
        super().__init__(sequence, TABLE_RNA_CODON_TO_AMINOACID)

    def get_sequence_map(self) -> dict:
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


    def dna_sequence(self) -> str:
        """
        Returns a string which represents the corresponding DNA sequence.
        """
        dna_sequence: str = self.sequence.replace("U", "T")
        return dna_sequence

    def to_dna(self):
        """
        Returns a DNA() object instance.
        """
        dna_sequence: str = self.dna_string()
        return DNA(dna_sequence)


class OpenReadFrame:
    def __init__(self, sequence: str, length_threshold: int = 0) -> None:
        """
        This class handles ORF related logic.\n
        This is not intended to be used without first creating a DNA instance.
        Whenever this class is instantiated, it will call the _get_orfs() method and store it as a local variable.
        By default, it will order the ORFs by size.
        """
        self.length_threshold = length_threshold
        self.fasta_map: list = self._get_fasta_sequence_map(sequence)
        self.orf_map: list = []
        self.update_orf_map()
        self.largest_fame: dict = {
        'nt_length': int(len(self.orf_map[0]['sequence'])),
        'aa_length': int(len(self.orf_map[0]['sequence']) / CODON_SIZE)
    }


    def _get_fasta_sequence_map(self, sequence: str) -> list:
        parser = FastaParserDNA()
        return parser.get_sequence_map(sequence)

    def extract_data_from_fasta_map(self,
            sequence: str,
            label: str,
            n: int,
            frame_reference: str,
            sequence_size: int) -> None:

        frame_sequence: str = ""
        
        start_codon: int = 0
        stop_codon: int = 0
        for i in range(sequence_size)[n::CODON_SIZE]:
            codon: str = sequence[i:i + CODON_SIZE]
            if frame_sequence == "":
                if codon == "ATG":
                    codon_start = self.find_start_codon(i, frame_reference, sequence_size)
                    frame_sequence += codon
            else:
                frame_sequence += codon
                if codon in STOP_CODON_DNA and frame_sequence != "":
                    if len(frame_sequence) > self.length_threshold:
                        values = {
                            'codon_start': codon_start,
                            'codon_stop': self.find_stop_codon(
                                i, frame_reference, sequence_size
                                ),
                            'count_aa': int( (len(frame_sequence) / CODON_SIZE) - 1),
                            'count_nt': int( (len(frame_sequence)) - 3),
                            'label': label,
                            'sequence': frame_sequence,
                            'reference': frame_reference
                            }

                        self.orf_map.append(values)
                        
                    frame_sequence = ""
                    # The frame should reset after a stop codon regardless of the size
                    # Only when the stop codon is found, then  the frame string will be appeneded to the list

    def find_start_codon(self,
        position: int,
        frame_reference: str,
        seq_size: int) -> int:
        if frame_reference in ["+", "++", "+++"]:
            position += 1
            return position
        elif frame_reference in ["-", "--", "---"]:
            position = abs(position - seq_size)
            return position

    def find_stop_codon(self,
        position: int,
        frame_reference: str,
        seq_size: int) -> int:
        if frame_reference in ["+", "++", "+++"]:
            position += 2
            return position
        elif frame_reference in ["-", "--", "---"]:
            position = abs((position + 2) - seq_size)
            return position


    def update_orf_map(self) -> dict:
        """        
        :length_threshold: int | minimal length of any given aminoacid sequence to be translated. Values inferior to this will be disqualified.\n

        This function will return a list with key-value pair relationship:\n
        'count_aa' is how many aminoacids are in the frame
        'count_nt' is the length of the frame
        'label' is the label of FASTA sequence determined by >
        'sequence' is the orf peptide sequence
        'reference' The orf_frame_reference:\n
            "+"   is the 1st possible read of the frame 5'->3'
            "++"  is the 2nd possible read of the frame 5'->3'
            "+++" is the 3rd possible read of the frame 5'->3'
            "-"   is the 1st possible read of the frame 3'->5'
            "--"  is the 2rd possible read of the frame 3'->5'
            "---" is the 3rd possible read of the frame 3'->5'


        """
        for item in self.fasta_map:
            dna = DNA(item["sequence"])
            label = item["label"]
            
            for i in range(CODON_SIZE*2):
                orf_frame_reference = ["+", "++", "+++", "-", "--", "---"]

                if i < CODON_SIZE:
                    self.extract_data_from_fasta_map(
                        sequence=dna.sequence,
                        label=label,
                        n=i,
                        frame_reference=orf_frame_reference[i],
                        sequence_size = dna.sequence_size
                    )

                elif i >= CODON_SIZE:
                    self.extract_data_from_fasta_map(
                        sequence=dna.template_strand(),
                        label=label,
                        n=(i-CODON_SIZE),
                        frame_reference=orf_frame_reference[i],
                        sequence_size = dna.sequence_size
                    )

        map_sorted = sorted(self.orf_map, key=lambda x: len(x['sequence']), reverse=True)
        
        return map_sorted
