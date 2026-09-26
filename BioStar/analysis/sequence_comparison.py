from BioStar.data.biochemistry import CODON_SIZE
from BioStar.engine.biochemistry import BiochemistryData
from BioStar.engine.default_data import get_default_biochemistry


class CompareNucleotideSequence:
    """Compares a reference DNA sequence against one or more sequences."""

    def __init__(
        self,
        original_sequence: str = "",
        compared_sequence: str | list[str] = "",
        data: BiochemistryData | None = None,
    ) -> None:
        self.original_sequence = original_sequence.upper()
        if isinstance(compared_sequence, str):
            self.compared_sequences = [compared_sequence.upper()]
        else:
            self.compared_sequences = [sequence.upper() for sequence in compared_sequence]
        self.data = data or get_default_biochemistry()
        self.classify_mutation = ClassifyNucleotideSequenceMutation(self.data)

    def compare(self, show_only_mutations: bool = True) -> list[dict[str, object]]:
        reference_codons = self._get_codons(self.original_sequence)
        results: list[dict[str, object]] = []

        for sequence in self.compared_sequences:
            test_codons = self._get_codons(sequence)
            self.classify_mutation.reset()

            for index, reference_codon in enumerate(reference_codons):
                if index >= len(test_codons):
                    break
                mutations = self.classify_mutation.classify(
                    reference_codon,
                    test_codons[index],
                )
                for mutation in mutations:
                    if show_only_mutations and mutation["mutation_type"] == "no_mutation":
                        continue
                    results.append(mutation)

        return results

    def _get_codons(self, sequence: str) -> list[str]:
        return [
            sequence[index:index + CODON_SIZE]
            for index in range(0, len(sequence) - (CODON_SIZE - 1), CODON_SIZE)
        ]


class ClassifyNucleotideSequenceMutation:
    """Classifies nucleotide-level mutations within DNA codons."""

    def __init__(self, data: BiochemistryData) -> None:
        self.codon_table = data.dna_codons
        self.stop_codons = data.dna_stop_codons
        self.reset()

    def reset(self) -> None:
        self.stop_codon_found = False
        self.codon_position = 1
        self.nucleotide_absolute_position = 1

    def classify(
        self,
        reference_codon: str,
        test_codon: str,
    ) -> list[dict[str, object]]:
        reference_aminoacid = self.codon_table[reference_codon]
        test_aminoacid = self.codon_table[test_codon]
        changed = reference_aminoacid != test_aminoacid

        if test_codon in self.stop_codons and changed:
            mutation_name = "nonsense"
            self.stop_codon_found = True
        elif changed:
            mutation_name = "missense"
        else:
            mutation_name = "no_mutation"

        return self.mutation_per_nucleotides_in_codon(
            mutation_name,
            reference_codon,
            test_codon,
            reference_aminoacid,
            test_aminoacid,
            changed,
        )

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
        changed_aminoacid: bool,
    ) -> dict[str, object]:
        return {
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
            "stop_codon_found_flag": self.stop_codon_found,
        }

    def mutation_per_nucleotides_in_codon(
        self,
        mutation_name: str,
        reference_codon: str,
        test_codon: str,
        reference_aminoacid: str,
        test_aminoacid: str,
        changed_aminoacid_flag: bool,
    ) -> list[dict[str, object]]:
        mutation_list: list[dict[str, object]] = []

        for index in range(CODON_SIZE):
            nucleotide_mutation = mutation_name
            if test_codon[index] == reference_codon[index]:
                nucleotide_mutation = "no_mutation"

            mutation_list.append(
                self.get_mutation_map(
                    new_nucleotide=test_codon[index],
                    old_nucleotide=reference_codon[index],
                    mutation=nucleotide_mutation,
                    nucleotide_absolute_position=self.nucleotide_absolute_position,
                    nucleotide_relative_position=index + 1,
                    new_codon=test_codon,
                    old_codon=reference_codon,
                    codon_position=self.codon_position,
                    changed_aminoacid=changed_aminoacid_flag,
                    new_aminoacid=test_aminoacid,
                    old_aminoacid=reference_aminoacid,
                )
            )
            self.nucleotide_absolute_position += 1

        self.codon_position += 1
        return mutation_list
