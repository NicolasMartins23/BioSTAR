from re import sub
from BioStar.DataBiochemistry import (AMINOACIDS, AMINOACID_TABLE,
    AMINOACIDS_AROMATIC, AMINOACIDS_NONPOLAR, AMINOACIDS_POLAR)

class Protein:
    def __init__(self, sequence: str = "") -> None:
        self.sequence: str = self._fasta_sequence(sequence)
        self.count: dict = self._get_count()
        self.sequence_size = self.count["total"]

    def _fasta_sequence(self, sequence: str = "") -> str:
        fasta_sequence = sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
        return fasta_sequence

    def _get_count(self) -> dict:
        count_by_aminoacid: dict = {aa: 0 for aa in AMINOACIDS}

        count_total: int = 0
        count_aromatic: int = 0
        count_nonpolar: int = 0
        count_polar: int = 0
        count_polar_negative: int = 0
        count_polar_neutral: int = 0
        count_polar_positive: int = 0

        for aa in self.sequence:
            if aa in AMINOACIDS_AROMATIC:
                count_by_aminoacid[aa] += 1
                count_aromatic += 1
                count_total += 1
            elif aa in AMINOACIDS_NONPOLAR:
                count_by_aminoacid[aa] += 1
                count_nonpolar += 1
                count_total += 1
            elif aa in AMINOACIDS_POLAR:
                count_by_aminoacid[aa] += 1
                count_polar += 1
                if aa in ["D", "E"]:
                    # Aspartic Acid and Glumatic Acid; respectively
                    count_polar_negative += 1
                elif aa in ["C", "N", "Q", "S", "T", "Y"]:
                    # Cystein, Asparagine, Glutamine, Serine, Threonine, Tyrosine
                    count_polar_neutral += 1
                elif aa in ["H", "K", "R"]:
                    count_polar_positive += 1
                    # Histidine, Lysine, Arginine; respectively

                count_total += 1

        
        aminoacid_count: dict = {
            'by_aminoacid': count_by_aminoacid,
            'aromatic': count_aromatic,
            'nonpolar': count_nonpolar,
            'polar': count_polar,
            'polar_negative': count_polar_negative,
            'polar_neutral': count_polar_neutral,
            'polar_positive': count_polar_positive,
            'total': count_total,
        }

        return aminoacid_count

    def get_aromacity(self) -> float:
        aromacity: float = self.count['polar'] / self.count['total']
        return round((aromacity * 100), 1)
