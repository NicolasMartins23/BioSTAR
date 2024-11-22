from re import sub
from BioStar.DataBiochemistry import (AMINOACIDS, AMINOACID_TABLE,
    AMINOACIDS_AROMATIC, AMINOACIDS_NONPOLAR,
    AMINOACIDS_POLAR, AMINOACIDS_POSITIVE, AMINOACIDS_NEGATIVE,
    C_TERM_PKA, N_TERM_PKA, WATER_MASS)

class Protein:
    def __init__(self, sequence: str = "") -> None:
        """
        Summary:
            Initializes the Protein object with a sequence and computes amino acid counts.

        Parameters:
            sequence (str): The protein sequence in FASTA format (default is an empty string).

        Returns:
            None
        """
        self.sequence: str = self._fasta_sequence(sequence)
        self.count: dict = self._update_count()
        self.sequence_size = self.count["total"]

    def _fasta_sequence(self, sequence: str = "") -> str:
        """
        Summary:
            Converts the input sequence into a valid FASTA format, removing non-amino acid characters.

        Parameters:
            sequence (str): The input protein sequence.

        Returns:
            str: The cleaned and uppercased sequence in FASTA format.
        """
        fasta_sequence = sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
        return fasta_sequence

    def _update_count(self) -> dict:
        """
        Summary:
            Updates and returns a dictionary with the counts of each type of amino acid in the protein sequence.

        Parameters:
            self

        Returns:
            dict: A dictionary containing counts for different categories of amino acids (e.g., aromatic, nonpolar, polar).
        """
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
                    count_polar_negative += 1
                elif aa in ["C", "N", "Q", "S", "T", "Y"]:
                    count_polar_neutral += 1
                elif aa in ["H", "K", "R"]:
                    count_polar_positive += 1
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
        """
        Summary:
            Calculates the aromacity of the protein, defined as the fraction of aromatic residues.

        Parameters:
            multiply_by (float): A multiplier to scale the result (default is 1.0).

        Returns:
            float: The aromacity of the protein, scaled by the multiplier.
        """
        aromacity: float = self.count['aromatic'] / self.count['total']
        return round((aromacity * multiply_by), 1)

    def charge_at_pH(self, pH: float) -> float:
        """
        Summary:
            Calculates the charge of the protein at a given pH based on its amino acid composition.

        Parameters:
            pH (float): The pH at which the protein charge is to be calculated.

        Returns:
            float: The net charge of the protein at the given pH, rounded to two decimal places.
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
        """
        Summary:
            Calculates the composition ratio of amino acids in the protein sequence.

        Parameters:
            multiply_by (float): A multiplier to scale the composition ratios (default is 1.0).

        Returns:
            dict: A dictionary of amino acid composition ratios, scaled by the multiplier.
        """
        total = self.count["total"]
        composition_summary = {aa: round((count / total) * multiply_by, 2) for aa, count in self.count['by_aminoacid'].items()}
        return composition_summary

    def extinction_coefficient(self) -> dict:
        """
        Summary:
            Calculates the extinction coefficient of the protein at 280 nm.
            Based on the number of Trp, Tyr, and Cys residues.

        Parameters:
            self

        Returns:
            dict: The extinction coefficient (in M^-1 cm^-1).
                cys_cystines assumes all cysteine residues form cystines.
                cys_reduced assumes all cysteine residues are reduced.
        """
        c_count = self.sequence.count("C")  # Cysteine
        w_count = self.sequence.count("W")  # Tryptophan
        y_count = self.sequence.count("Y")  # Tyrosine

        # Assuming all Cys form disulfide bonds
        disulfide_bonds = c_count // 2

        # Extinction coefficients (in M^-1 cm^-1)
        c_coeff = 125
        w_coeff = 5500
        y_coeff = 1490

        # Total extinction coefficient
        coeff_no_c = (
            (w_count * w_coeff) +
            (y_count * y_coeff)
        )

        coeff_with_c = coeff_no_c + (disulfide_bonds * c_coeff)

        coefficients: dict = {
            "cys_cystines": round(coeff_with_c, 2),
            "cys_reduced": round(coeff_no_c, 2),
        }

        return coefficients

    def hydrophobic_index(self) -> float:
        """
        Summary:
            Calculates the hydrophobic index of the protein based on the hydrophobicity of each amino acid.

        Parameters:
            self

        Returns:
            float: The average hydrophobicity of the protein, rounded to two decimal places.
        """
        h_index = sum(AMINOACID_TABLE[aa]['hydrophobicity'] for aa in self.sequence)
        return round(h_index / len(self.sequence), 2)

    def isoelectric_point(self) -> float:
        """
        Summary:
            Calculates the isoelectric point (pI) of the protein, the pH at which the net charge of the protein is zero.

        Parameters:
            self

        Returns:
            float: The estimated isoelectric point (pI) of the protein, rounded to two decimal places.
        """
        low_pH = 0.0
        high_pH = 14.0
        tolerance = 0.01  # Acceptable deviation for "zero charge"

        while high_pH - low_pH > 0.01:  # Precision up to 2 decimal places
            mid_pH = (low_pH + high_pH) / 2.0
            net_charge = self.charge_at_pH(mid_pH)

            if abs(net_charge) < tolerance:  # Net charge close enough to zero
                return round(mid_pH, 2)
            elif net_charge > 0:  # Protein is positively charged
                low_pH = mid_pH
            else:  # Protein is negatively charged
                high_pH = mid_pH

        # Return the midpoint as the best estimate of the pI
        return round((low_pH + high_pH) / 2.0, 2)

    def molecular_weight(self) -> float:
        """
        Summary:
            Calculates the molecular weight of the protein by summing the weights of its amino acids.

        Parameters:
            self

        Returns:
            float: The molecular weight of the protein, with the mass of water subtracted.
        """
        molecular_weight = sum(AMINOACID_TABLE[aa]['weight'] for aa in self.sequence)

        # Subtract water from the amount of peptide bonds
        molecular_weight -= (self.count["total"] - 1) * WATER_MASS

        return round(molecular_weight, 2)

    def secondary_structure_propensity(self) -> dict:
        """
        Summary:
            Calculates the propensity for secondary structures (alpha helix, beta sheet, coil) in the protein.

        Parameters:
            self

        Returns:
            dict: A dictionary containing the percentage propensity for alpha helix, beta sheet, and coil.
        """
        alpha_helix_propensity = sum(AMINOACID_TABLE[aa]['alpha_helix'] for aa in self.sequence) / len(self.sequence)
        beta_sheet_propensity = sum(AMINOACID_TABLE[aa]['beta_sheet'] for aa in self.sequence) / len(self.sequence)
        coil_propensity = 1 - (alpha_helix_propensity + beta_sheet_propensity)  # Simplified assumption
        return {
            'alpha_helix': round(alpha_helix_propensity * 100, 2),
            'beta_sheet': round(beta_sheet_propensity * 100, 2),
            'coil': round(coil_propensity * 100, 2)
        }
