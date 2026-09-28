from __future__ import annotations

import gzip
import sqlite3
from functools import lru_cache
from importlib.resources import files

from BioStar.engine.biochemistry import AminoAcidData, BiochemistryData


@lru_cache(maxsize=1)
def get_default_biochemistry() -> BiochemistryData:
    database_bytes = gzip.decompress(
        files("BioStar.data").joinpath("biostar.sqlite3.gz").read_bytes()
    )
    connection = sqlite3.connect(":memory:")
    try:
        connection.deserialize(database_bytes)
        return BiochemistryData(
            amino_acids=_load_amino_acids(connection),
            aromatic=_load_class(connection, "aromatic"),
            nonpolar=_load_class(connection, "nonpolar"),
            polar=_load_class(connection, "polar"),
            positive=_load_class(connection, "positive"),
            negative=_load_class(connection, "negative"),
            dna_codons=_load_codons(connection, "DNA"),
            rna_codons=_load_codons(connection, "RNA"),
            dna_stop_codons=_load_stop_codons(connection, "DNA"),
            rna_stop_codons=_load_stop_codons(connection, "RNA"),
            dna_start_codon=_load_start_codon(connection, "DNA"),
            rna_start_codon=_load_start_codon(connection, "RNA"),
            water_mass=_load_constant(connection, "water_mass"),
            n_term_pka=_load_constant(connection, "n_term_pka"),
            c_term_pka=_load_constant(connection, "c_term_pka"),
        )
    finally:
        connection.close()


def _load_amino_acids(connection: sqlite3.Connection) -> dict[str, AminoAcidData]:
    rows = connection.execute(
        """
        SELECT symbol, molecular_weight, hydrophobicity, alpha_helix,
               beta_sheet, pka, pkb, pkr
        FROM amino_acids
        WHERE symbol != '*'
        ORDER BY symbol
        """
    )
    return {
        row[0]: AminoAcidData(
            symbol=row[0],
            molecular_weight=float(row[1]),
            hydrophobicity=float(row[2]),
            alpha_helix=float(row[3]),
            beta_sheet=float(row[4]),
            pka=float(row[5]) if row[5] is not None else None,
            pkb=float(row[6]) if row[6] is not None else None,
            pkr=float(row[7]) if row[7] is not None else None,
        )
        for row in rows
    }


def _load_class(connection: sqlite3.Connection, class_name: str) -> frozenset[str]:
    rows = connection.execute(
        """
        SELECT amino_acid_symbol
        FROM amino_acid_class_members
        WHERE class_name = ?
        """,
        (class_name,),
    )
    return frozenset(row[0] for row in rows)


def _load_codons(connection: sqlite3.Connection, molecule_type: str) -> dict[str, str]:
    rows = connection.execute(
        """
        SELECT sequence, amino_acid_symbol
        FROM codons
        WHERE genetic_code_id = 1 AND molecule_type = ?
        """,
        (molecule_type,),
    )
    return {row[0]: row[1] for row in rows}


def _load_stop_codons(connection: sqlite3.Connection, molecule_type: str) -> frozenset[str]:
    rows = connection.execute(
        """
        SELECT sequence
        FROM codons
        WHERE genetic_code_id = 1 AND molecule_type = ? AND is_stop = 1
        """,
        (molecule_type,),
    )
    return frozenset(row[0] for row in rows)


def _load_start_codon(connection: sqlite3.Connection, molecule_type: str) -> str:
    row = connection.execute(
        """
        SELECT sequence
        FROM codons
        WHERE genetic_code_id = 1 AND molecule_type = ? AND is_start = 1
        """,
        (molecule_type,),
    ).fetchone()
    if row is None:
        raise RuntimeError(f"No start codon found for {molecule_type}.")
    return str(row[0])


def _load_constant(connection: sqlite3.Connection, key: str) -> float:
    row = connection.execute(
        "SELECT value FROM reference_constants WHERE key = ?",
        (key,),
    ).fetchone()
    if row is None:
        raise RuntimeError(f"Missing biochemical reference constant: {key}.")
    return float(row[0])
