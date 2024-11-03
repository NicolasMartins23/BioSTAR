from ..DataBiochemistry import *

class CompareNucleotideSequence:
    """
    This class for compares two or more nucleotide sequences
    and identifies mutations that may occur in from one sequence to another.
    """
    def __init__(self,
                original_sequence: list[str] = [""],
                compared_sequence: list[str] = [""]) -> None:
        self.original_sequence: list[str] = original_sequence
        self.compared_sequence: list[str] = compared_sequence
        self.classify_mutation: ClassifyNucleotideSequenceMutation = ClassifyNucleotideSequenceMutation()

    def compare(self, show_only_mutations: bool = True):
        """
        This is the main method of the CompareNucleotide class.\n
        It will compare a reference sequence to a test sequence list\n
        The algorithm requires that the test sequence list is in fasta format\n
        It will return a dictonary containing the differences between each aminoacid
        in the test sequence list with the following structure:\n
        result = [
            {\n
                "nucleotide": "A" "C" "G" "T" or "U",\n
                "mutation_type": "no_mutation", "silent", "missense", "nonsense", "nontranslated"\n
                "changed_aminoacid_flag": True or False
            },\n
            (...) for every nucleotide
        ]
        no_mutation = The nucleotide is the same.\n
        silent = The nucleotide changed, but it translates the same aminoacid.\n
        missense_dg = The nucleotide changed and the codon translates a different aminoacid of a different group (nonpolar or polar).\n
        missense_sg = The nucleotide changed and the codon translates a different aminoacid of the samge group (nonpolar or polar).\n
        nonsense = The nucleotide changed and the codon translates a stop codon .\n
        not_translated = The nucleotide is not translated. It is the mutation_type for every nucleotide after a nonsense translation happens.
        """

        reference_codon_values = self._get_reference_codon_values()
        test_codon_values = self._get_test_codon_values()

        return self._search_for_mutations(
            reference_codon_values,
            test_codon_values,
            show_only_mutations)


    def _get_test_codon_values(self) -> list[str]:
        values: list[str] = []
        for sequence in [self.compared_sequence]:
            codons: list = []
            for i in range(0, len(sequence), CODON_SIZE):
                codons.append(sequence[i:i + CODON_SIZE])
            values.append(codons)

        return values

    def _get_reference_codon_values(self) -> list[str]:
        values: list[str] = []
        for i in range(0, len(self.original_sequence), CODON_SIZE):
            values.append(self.original_sequence[i:i + CODON_SIZE])

        return values

    def _search_for_mutations(self, referece_sequence_list, compared_sequence, show_only_mutations) -> list[dict]:
        """
       result: [
            {
                "nucleotide': "A", "C", "G", "T", or "U"
                "mutation_type": "",
                "changed_aminoacid_flag": True or False
            }
        ]
        """
        seq_size: int = len(referece_sequence_list)
        result: list[dict] = []

        for test_sequence in compared_sequence:
            for i in range(seq_size):
                reference_codon: str = referece_sequence_list[i]
                test_codon: str = test_sequence[i]

                mutations = self.classify_mutation.classify(
                    reference_codon, test_codon)
                for mutation in mutations:
                    if show_only_mutations:
                        if mutation['mutation_type'] != "no_mutation":
                            result.append(mutation)
                    else:
                        result.append(mutation)
        return result


class ClassifyNucleotideSequenceMutation:
    """
    This class classifies the different mutation types when a mutation is identified.
    It analyses the mutations for every codon
    """
    def __init__(self) -> None:
        self.CODON_TABLE_REF: list = TABLE_DNA_CODON_TO_AMINOACID
        self.STOP_CODON_REF: list = STOP_CODON_DNA
        self.stop_codon_found: bool = False
        self.__codon_position: int = 1
        self.__nucleotide_absolute_position: int = 1

    def classify(
        self,
        reference_codon: str,
        test_codon: str) -> list:
        """
        This function will return a mutation list with mutation details for
        every nucleotide in the codon.
        """

        reference_aminoacid = self.CODON_TABLE_REF[reference_codon]
        test_aminoacid = self.CODON_TABLE_REF[test_codon]

        aminoacid_has_changed: bool = (reference_aminoacid != test_aminoacid)

        mutation_results: list = []

        if aminoacid_has_changed:

            if test_codon in self.STOP_CODON_REF:
                nonsense_codon_map = self._mutation_nonsense(
                    reference_codon, test_codon,
                    reference_aminoacid, test_aminoacid
                    )
                for nucleotide in nonsense_codon_map:
                    mutation_results.append(nucleotide)

                self.stop_codon_found = True

            else:
                missense_codon_map = self._mutation_missense(
                    reference_codon, test_codon,
                    reference_aminoacid, test_aminoacid
                    )
                for nucleotide in missense_codon_map:
                    mutation_results.append(nucleotide)

        else:
            missense_codon_map = self._mutation_silent(
                reference_codon, test_codon,
                reference_aminoacid, test_aminoacid
                )
            for nucleotide in missense_codon_map:
                mutation_results.append(nucleotide)

        return mutation_results

    def get_mutation_map(
            self,
            new_nucleotide: str,
            old_nucleotide: str,
            mutation: str,
            new_codon: str,
            old_codon: str,
            new_aminoacid: str,
            old_aminoacid: str,
            codon_position: int,
            nucleotide_absolute_position: int,
            nucleotide_relative_position: int,
            changed_aminoacid: bool
            ) -> dict:
        nucleotide_mutation_map = {
            "mutation_type": mutation,
            "nucleotide_details": {
                "new_nucleotide": new_nucleotide,
                "old_nucleotide": old_nucleotide,
                "nucleotide_absolute_position": nucleotide_absolute_position,
                "nucleotide_relative_position": nucleotide_relative_position,

            },
            "codon_details": {
                "new_codon": new_codon,
                "old_codon": old_codon,
                "codon_position": codon_position,
            },
            "aminoacid_details": {
                "new_aminoacid": new_aminoacid,
                "old_aminoacid": old_aminoacid,
                "changed_aminoacid_flag": changed_aminoacid,
            },
            
            "stop_codon_found_flag": self.stop_codon_found
        }

        return nucleotide_mutation_map

    def mutation_per_nucleotides_in_codon(
            self,
            mutation_name: str,
            reference_codon: str,
            test_codon: str,
            reference_aminoacid: str,
            test_aminoacid: str,
            changed_aminoacid_flag: bool,
        ) -> list[dict]:

        mutation_list = []


        for i in range(CODON_SIZE):
            if test_codon[i] == reference_codon[i]:
                mutation_list.append(
                    self.get_mutation_map(
                        new_nucleotide=test_codon[i],
                        old_nucleotide=reference_codon[i],
                        mutation="no_mutation",
                        nucleotide_absolute_position=self.__nucleotide_absolute_position,
                        nucleotide_relative_position=i+1,
                        new_codon=test_codon,
                        old_codon=reference_codon,
                        codon_position=self.__codon_position,
                        changed_aminoacid=changed_aminoacid_flag,
                        new_aminoacid=test_aminoacid,
                        old_aminoacid=reference_aminoacid,
                    )
                )
            else:
                mutation_list.append(
                    self.get_mutation_map(
                        new_nucleotide=test_codon[i],
                        old_nucleotide=reference_codon[i],
                        mutation=mutation_name,
                        nucleotide_absolute_position=self.__nucleotide_absolute_position,
                        nucleotide_relative_position=i+1,
                        new_codon=test_codon,
                        old_codon=reference_codon,
                        codon_position=self.__codon_position,
                        changed_aminoacid=changed_aminoacid_flag,
                        new_aminoacid=test_aminoacid,
                        old_aminoacid=reference_aminoacid,
                    )
                )

            self.__nucleotide_absolute_position += 1  # Updates the nucleotide counter
        
        self.__codon_position +=1  # Updates the nucleotide counter
        return mutation_list

    def _mutation_silent(
            self,
            reference_codon,
            test_codon,
            reference_aminoacid,
            test_aminoacid
        ) -> list[dict]:

        return self.mutation_per_nucleotides_in_codon(
            mutation_name="silent",
            reference_codon=reference_codon,
            test_codon=test_codon,
            reference_aminoacid=reference_aminoacid,
            test_aminoacid=test_aminoacid,
            changed_aminoacid_flag=True
        )

    def _mutation_nonsense(
            self,
            reference_codon,
            test_codon,
            reference_aminoacid,
            test_aminoacid
        ) -> list[dict]:

        return self.mutation_per_nucleotides_in_codon(
            mutation_name="nonsense",
            reference_codon=reference_codon,
            test_codon=test_codon,
            reference_aminoacid=reference_aminoacid,
            test_aminoacid=test_aminoacid,
            changed_aminoacid_flag=True
        )

    def _mutation_missense(
            self,
            reference_codon,
            test_codon,
            reference_aminoacid,
            test_aminoacid
        ) -> list[dict]:

        return self.mutation_per_nucleotides_in_codon(
            mutation_name="missense",
            reference_codon=reference_codon,
            test_codon=test_codon,
            reference_aminoacid=reference_aminoacid,
            test_aminoacid=test_aminoacid,
            changed_aminoacid_flag=True
        )
