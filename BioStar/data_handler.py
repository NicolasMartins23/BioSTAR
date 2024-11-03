class FastaParser:
    def __init__(self) -> None:
        self.sequence: list[str] = [""]

    def get_sequence_map(
            self,
            fasta_string: str,
            ) -> list:

        sequence_map: list = []
        lines: list = fasta_string.strip().split('\n')

        sequence_map_entry = None


        for line in lines:
            if line.startswith('>'):
                if sequence_map_entry:
                    sequence_map.append(sequence_map_entry)
                sequence_map_entry = {
                    "label": line[1:].strip(),
                    "sequence": ""
                }
            elif sequence_map_entry:
                sequence_map_entry["sequence"] += line.strip()
            else:
                # Handle the case where the first line does not start with '>'
                if sequence_map_entry is None:
                    sequence_map_entry = {
                        "label": f"MySeq",
                        "sequence": line.strip()
                    }
                else:
                    sequence_map_entry["sequence"] += line.strip()

        if sequence_map_entry:
            sequence_map.append(sequence_map_entry)

        return sequence_map


class FastaParserDNA(FastaParser):
    pass
