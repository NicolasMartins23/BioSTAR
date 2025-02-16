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
        Translates the nucleic acid sequence into a peptide sequence and returns it as a string.
    
        Parameters:
            show_stop_codon (bool): If True, includes the stop codon in the returned peptide sequence. 
                                    Defaults to False, which excludes the stop codon.
    
        Returns:
            str: The peptide sequence derived from the nucleic acid, translated based on the codon table. 
                 If show_stop_codon is False, the sequence excludes the stop codon.
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
        Converts the nucleic acid sequence into a Protein object.
    
        Returns:
            Protein: A Protein object containing the peptide sequence translated from the nucleic acid sequence.
        """
        peptide_sequence = self.get_peptide_sequence()
        return Protein(peptide_sequence)


class DNA(NucleicAcid):
    """
    Handles DNA-specific functions and inherits NucleicAcid functionality.
    It must be initialized with a valid nucleotide sequence as a parameter.
    """
    def __init__(self, sequence: str) -> None:
        """
        Initializes a DNA object with a nucleotide sequence and sets the codon table for DNA.

        Parameters:
            sequence (str): A string representing the DNA sequence.
        """
        super().__init__(sequence, TABLE_DNA_CODON_TO_AMINOACID)

    def get_sequence_map(self) -> dict:
        """
        Calculates the frequency of each nucleotide (A, C, G, T) in the DNA sequence.

        Returns:
            dict: A dictionary with counts for each nucleotide and the total count.
        """
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
        """
        Converts a given DNA sequence to FASTA format by removing invalid characters.

        Parameters:
            dna_seq (str): A DNA sequence string. Defaults to an empty string.

        Returns:
            str: A valid DNA sequence containing only A, C, G, and T.
        """
        fasta_sequence = re.sub(r'[^ACGT]', '', dna_seq.upper())
        return fasta_sequence

    def gc_content(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        """
        Calculates the GC content ratio of the DNA sequence.

        Parameters:
            multiply_by (float): A multiplier for the GC content value. Defaults to 1.0.

        Returns:
            float: The GC content ratio, rounded to two decimal places.
        """
        gc_count: int = self.sequence_map["C"] + self.sequence_map["G"]
        gc_content: float = (gc_count / self.sequence_size) * multiply_by

        return round(gc_content, decimal_places)

    def at_skew(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        """
        Calculates the skew of Adenine (A) and Thymine (T) content in the DNA sequence.

        Parameters:
            multiply_by (float): A multiplier for the skew value. Defaults to 1.0.

        Returns:
            float: The AT skew, where a higher value indicates an Adenine prevalence,
                   and a lower value indicates a Thymine prevalence, rounded to two decimal places.
        """
        at_skew: float = (self.sequence_map["A"] - self.sequence_map["T"]) / (self.sequence_map["A"] + self.sequence_map["T"])

        return round((at_skew * multiply_by), decimal_places)

    def gc_skew(self, multiply_by: float = 1.0, decimal_places: int = 4) -> float:
        """
        Calculates the skew of Guanine (G) and Cytosine (C) content in the DNA sequence.

        Parameters:
            multiply_by (float): A multiplier for the skew value. Defaults to 1.0.

        Returns:
            float: The GC skew, where a higher value indicates a Guanine prevalence,
                   and a lower value indicates a Cytosine prevalence, rounded to two decimal places.
        """
        gc_skew: float = (self.sequence_map["G"] - self.sequence_map["C"]) / (self.sequence_map["G"] + self.sequence_map["C"])
        return round((gc_skew * multiply_by), decimal_places)

    def template_strand(self, reverse_string: bool = True) -> str:
        """
        Generates the complementary strand of the DNA sequence (template strand).

        Parameters:
            reverse_string (bool): If True, returns the reversed complementary strand. Defaults to True.

        Returns:
            str: The template strand of the DNA sequence.
        """
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
        Converts the DNA sequence into the corresponding RNA sequence.

        Returns:
            str: A string representing the RNA sequence, with Thymine (T) replaced by Uracil (U).
        """
        rna_sequence: str = self.sequence.replace("T", "U")
        return rna_sequence

    def to_rna(self):
        """
        Converts the DNA sequence into an RNA object.

        Returns:
            RNA: An RNA object created from the converted RNA sequence.
        """
        rna_sequence: str = self.rna_sequence()
        return RNA(rna_sequence)

    def orf_map(self, length_threshold: int = 0) -> dict:
        """
        Identifies open reading frames (ORFs) in the DNA sequence.

        Parameters:
            length_threshold (int): Minimum length for an ORF to be included. Defaults to 0.

        Returns:
            dict: A dictionary representing the identified ORFs and their positions in the sequence.
        """
        orf = OpenReadFrame(sequence=self.sequence, length_threshold=length_threshold)
        return orf.orf_map


class RNA(NucleicAcid):
    """
    Handles RNA-specific functions and inherits NucleicAcid functionality.
    It must be initialized with a valid RNA nucleotide sequence as a parameter.
    """
    def __init__(self, sequence: str) -> None:
        """
        Initializes an RNA object with a nucleotide sequence and sets the codon table for RNA.

        Parameters:
            sequence (str): A string representing the RNA sequence.
        """
        super().__init__(sequence, TABLE_RNA_CODON_TO_AMINOACID)

    def get_sequence_map(self) -> dict:
        """
        Calculates the frequency of each nucleotide (A, C, G, U) in the RNA sequence.

        Returns:
            dict: A dictionary with counts for each nucleotide and the total count.
        """
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
        """
        Converts a given RNA sequence to FASTA format by removing invalid characters.

        Parameters:
            sequence (str): An RNA sequence string. Defaults to an empty string.

        Returns:
            str: A valid RNA sequence containing only A, C, G, and U.
        """
        fasta_sequence = re.sub(r'[^ACGU]', '', sequence.upper())
        return fasta_sequence

    def trim_on_stop_codon(self) -> str:
        """
        Trims the RNA sequence at the first stop codon.

        Returns:
            str: The RNA sequence truncated at the first stop codon. If no stop codon is found, 
                 the full sequence is returned.
        """
        trimmed_sequence: str = ""
        for i in range(0, self.sequence_size, CODON_SIZE):
            codon: str = self.sequence[i:i + CODON_SIZE]
            trimmed_sequence += codon
            if codon in STOP_CODON_RNA:
                return trimmed_sequence
        return trimmed_sequence

    def dna_sequence(self) -> str:
        """
        Converts the RNA sequence into the corresponding DNA sequence.

        Returns:
            str: A string representing the DNA sequence, with Uracil (U) replaced by Thymine (T).
        """
        dna_sequence: str = self.sequence.replace("U", "T")
        return dna_sequence

    def to_dna(self):
        """
        Converts the RNA sequence into a DNA object.

        Returns:
            DNA: A DNA object created from the converted DNA sequence.
        """
        dna_sequence: str = self.dna_sequence()
        return DNA(dna_sequence)


class OpenReadFrame:
    """
    Handles Open Reading Frame (ORF) detection and related logic.
    This class is designed to be used in conjunction with a DNA object.
    When instantiated, it initializes by parsing the provided sequence into FASTA format,
    calculating ORFs, and sorting them by size.
    """
    def __init__(self, sequence: str, length_threshold: int = 0) -> None:
        """
        Initializes an OpenReadFrame instance with a DNA sequence.

        Parameters:
            sequence (str): The nucleotide sequence to process.
            length_threshold (int): Minimum length of an ORF to include in results. Defaults to 0.
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
        """
        Converts a sequence into a FASTA-formatted sequence map.

        Parameters:
            sequence (str): The nucleotide sequence to convert.

        Returns:
            list: A list of parsed FASTA sequences.
        """
        parser = FastaParserDNA()
        return parser.get_sequence_map(sequence)

    def extract_data_from_fasta_map(
        self,
        sequence: str,
        label: str,
        n: int,
        frame_reference: str,
        sequence_size: int
    ) -> None:
        """
        Extracts ORF data from a single FASTA-formatted sequence and updates the ORF map.

        Parameters:
            sequence (str): The nucleotide sequence to analyze.
            label (str): The label associated with the FASTA sequence.
            n (int): Frame offset (start position).
            frame_reference (str): Frame reference indicator (e.g., "+", "-", "++").
            sequence_size (int): The length of the nucleotide sequence.
        """
        frame_sequence: str = ""
        codon_start: int = 0

        for i in range(sequence_size)[n::CODON_SIZE]:
            codon: str = sequence[i:i + CODON_SIZE]
            if frame_sequence == "":
                if codon == "ATG":  # Start codon
                    codon_start = self.find_start_codon(i, frame_reference, sequence_size)
                    frame_sequence += codon
            else:
                frame_sequence += codon
                if codon in STOP_CODON_DNA and frame_sequence != "":
                    if len(frame_sequence) > self.length_threshold:
                        values = {
                            'codon_start': codon_start,
                            'codon_stop': self.find_stop_codon(i, frame_reference, sequence_size),
                            'count_aa': int((len(frame_sequence) / CODON_SIZE) - 1),
                            'count_nt': int(len(frame_sequence) - 3),
                            'label': label,
                            'sequence': frame_sequence,
                            'reference': frame_reference
                        }
                        self.orf_map.append(values)
                    frame_sequence = ""  # Reset after a stop codon

    def find_start_codon(
        self,
        position: int,
        frame_reference: str,
        seq_size: int
    ) -> int:
        """
        Finds the start position of a codon in the specified reading frame.

        Parameters:
            position (int): Current position in the sequence.
            frame_reference (str): Frame reference (e.g., "+", "-", "++").
            seq_size (int): Total sequence size.

        Returns:
            int: Adjusted start position of the codon.
        """
        if frame_reference in ["+", "++", "+++"]:
            position += 1
            return position
        elif frame_reference in ["-", "--", "---"]:
            position = abs(position - seq_size)
            return position

    def find_stop_codon(
        self,
        position: int,
        frame_reference: str,
        seq_size: int
    ) -> int:
        """
        Finds the stop position of a codon in the specified reading frame.

        Parameters:
            position (int): Current position in the sequence.
            frame_reference (str): Frame reference (e.g., "+", "-", "++").
            seq_size (int): Total sequence size.

        Returns:
            int: Adjusted stop position of the codon.
        """
        if frame_reference in ["+", "++", "+++"]:
            position += 2
            return position
        elif frame_reference in ["-", "--", "---"]:
            position = abs((position + 2) - seq_size)
            return position

    def update_orf_map(self) -> dict:
        """
        Updates the ORF map by identifying potential ORFs in all reading frames.

        Returns:
            dict: A list of identified ORFs, sorted by sequence length.

        Key-value pairs in the ORF map:
            - 'count_aa': Number of amino acids in the ORF.
            - 'count_nt': Length of the ORF in nucleotides.
            - 'label': Label of the FASTA sequence (e.g., ">seq1").
            - 'sequence': ORF peptide sequence.
            - 'reference': Reading frame reference:
                "+"   - 1st frame (5'->3')
                "++"  - 2nd frame (5'->3')
                "+++" - 3rd frame (5'->3')
                "-"   - 1st frame (3'->5')
                "--"  - 2nd frame (3'->5')
                "---" - 3rd frame (3'->5')
        """
        for item in self.fasta_map:
            dna = DNA(item["sequence"])
            label = item["label"]

            for i in range(CODON_SIZE * 2):
                orf_frame_reference = ["+", "++", "+++", "-", "--", "---"]

                if i < CODON_SIZE:
                    self.extract_data_from_fasta_map(
                        sequence=dna.sequence,
                        label=label,
                        n=i,
                        frame_reference=orf_frame_reference[i],
                        sequence_size=dna.sequence_size
                    )
                elif i >= CODON_SIZE:
                    self.extract_data_from_fasta_map(
                        sequence=dna.template_strand(),
                        label=label,
                        n=(i - CODON_SIZE),
                        frame_reference=orf_frame_reference[i],
                        sequence_size=dna.sequence_size
                    )

        map_sorted = sorted(self.orf_map, key=lambda x: len(x['sequence']), reverse=True)
        return map_sorted

