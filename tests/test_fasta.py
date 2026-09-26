from BioStar.io import FastaParser


def test_fasta_parser() -> None:
    parser = FastaParser()
    result = parser.get_sequence_map(">seq1\nATGC\n>seq2\nGG")
    assert result == [
        {"label": "seq1", "sequence": "ATGC"},
        {"label": "seq2", "sequence": "GG"},
    ]
