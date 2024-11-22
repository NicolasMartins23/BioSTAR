from re import sub
from BioStar.DataBiochemistry import (AMINOACIDS, AMINOACID_TABLE,
    AMINOACIDS_AROMATIC, AMINOACIDS_NONPOLAR,
    AMINOACIDS_POLAR, AMINOACIDS_POSITIVE, AMINOACIDS_NEGATIVE,
    C_TERM_PKA, N_TERM_PKA, WATER_MASS)

class Protein:
    def __init__(self, sequence: str = "") -> None:
        self.sequence: str = self._fasta_sequence(sequence)
        self.count: dict = self._update_count()
        self.sequence_size = self.count["total"]

    def _fasta_sequence(self, sequence: str = "") -> str:
        fasta_sequence = sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
        return fasta_sequence

    def _update_count(self) -> dict:
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

    def aromacity(self, multiply_by: float = 1.0) -> float:
        aromacity: float = self.count['aromatic'] / self.count['total']
        return round((aromacity * multiply_by), 1)

    def charge_at_ph(self, pH: float) -> float:
        """
        This function takes the pH as a parameter.
        Returns the floating point value of the estimated charge of a protein.
        """
        normalized_pH: float = round(pH, 3)
        positive_charge: float = 0.0
        negative_charge: float = 0.0

        # Charged residues and their pKa values
        positive_residues = AMINOACIDS_POSITIVE # Lys, Arg, His
        negative_residues = AMINOACIDS_NEGATIVE  # Asp, Glu, Cys, Tyr

        # Positive residues
        for aa in positive_residues:
            if aa in self.count['by_aminoacid']:
                count = self.count['by_aminoacid'][aa]
                pKa = AMINOACID_TABLE[aa]['pKr']
                if pKa is not None:  # Ensure pKa is defined
                    partial_charge = 1.0 / (1.0 + 10 ** (normalized_pH - pKa))
                    positive_charge += count * partial_charge

        # Negative residues
        for aa in negative_residues:
            if aa in self.count['by_aminoacid']:
                count = self.count['by_aminoacid'][aa]
                pKa = AMINOACID_TABLE[aa]['pKr']
                if pKa is not None:  # Ensure pKa is defined
                    partial_charge = 1.0 / (1.0 + 10 ** (pKa - normalized_pH))
                    negative_charge += count * partial_charge

        # N-terminal contribution
        partial_charge = 1.0 / (1.0 + 10 ** (normalized_pH - N_TERM_PKA))
        positive_charge += partial_charge

        # C-terminal contribution
        partial_charge = 1.0 / (1.0 + 10 ** (C_TERM_PKA - normalized_pH))
        negative_charge += partial_charge

        # Net charge
        return round(positive_charge - negative_charge, 2)


    def composition_ratio(self, multiply_by: float = 1.0) -> dict:
        '''
        Method which returns the ratio composition of all aminoacids.
        By default the sum adds up to 1 but you can change the multiplier.
        '''
        total = self.count["total"]
        composition_summary = {aa: round((count / total) * multiply_by, 2) for aa, count in self.count['by_aminoacid'].items()}
        return composition_summary

    def hydrophobic_index(self) -> float:
        '''
        Hydrophic index calculated by the sum of the hydrophobicity
        of every aminoacid present in the protein
        '''
        h_index = sum(AMINOACID_TABLE[aa]['hydrophobicity'] for aa in self.sequence)
        return round(hydrophobicity / len(self.sequence), 2)

    def molecular_weight(self) -> float:
        '''
        Sum the weight of every aminoacid present in the protein
        '''
        
        molecular_weight = sum(AMINOACID_TABLE[aa]['weight'] for aa in self.sequence)

        # Subtract water from the amount of peptide bonds
        molecular_weight -= (self.count["total"] - 1) * WATER_MASS

        return round(molecular_weight, 2)

    def secondary_structure_propensity(self) -> dict:
        '''
        Function returns a dictionary of alpha_helix, beta_sheet or coil propensity

        Propensity for secondary structures represents an intrinsic property of amino acid,
        and it is used for generating new algorithms and prediction methods.
        '''

        alpha_helix_propensity = sum(AMINOACID_TABLE[aa]['alpha_helix'] for aa in self.sequence) / len(self.sequence)
        beta_sheet_propensity = sum(AMINOACID_TABLE[aa]['beta_sheet'] for aa in self.sequence) / len(self.sequence)
        coil_propensity = 1 - (alpha_helix_propensity + beta_sheet_propensity)  # Simplified assumption
        return {
            'alpha_helix': round(alpha_helix_propensity * 100, 2),
            'beta_sheet': round(beta_sheet_propensity * 100, 2),
            'coil': round(coil_propensity * 100, 2)
        }
