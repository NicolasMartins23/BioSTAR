from __future__ import annotations

import gzip
import os
import sqlite3
from functools import lru_cache
from importlib.resources import files
from tempfile import NamedTemporaryFile

from BioStar.engine.biochemistry import AminoAcidData


class ReferenceData:
    """Provides read-only access to bundled biochemical reference data."""

    def __init__(self) -> None:
        database_bytes = gzip.decompress(
            files("BioStar.data").joinpath("biostar.sqlite3.gz").read_bytes()
        )
        with NamedTemporaryFile(suffix=".sqlite3", delete=False) as database_file:
            database_file.write(database_bytes)
            self._database_path: str = database_file.name
        self._connection: sqlite3.Connection = sqlite3.connect(
            self._database_path
        )

    def amino_acid(self, symbol: str) -> AminoAcidData:
        row = self._connection.execute(
            """
            SELECT symbol, molecular_weight, hydrophobicity, alpha_helix,
                   beta_sheet, pka, pkb, pkr
            FROM amino_acids
            WHERE symbol = ?
            """,
            (symbol,),
        ).fetchone()
        if row is None:
            raise KeyError(f"Unknown amino acid: {symbol}.")
        return AminoAcidData(
            symbol=str(row[0]),
            molecular_weight=float(row[1]),
            hydrophobicity=float(row[2]),
            alpha_helix=float(row[3]),
            beta_sheet=float(row[4]),
            pka=float(row[5]) if row[5] is not None else None,
            pkb=float(row[6]) if row[6] is not None else None,
            pkr=float(row[7]) if row[7] is not None else None,
        )

    def amino_acid_symbols(self) -> tuple[str, ...]:
        rows = self._connection.execute(
            "SELECT symbol FROM amino_acids WHERE symbol != '*' ORDER BY symbol"
        )
        return tuple(str(row[0]) for row in rows)

    def amino_acid_class(self, class_name: str) -> frozenset[str]:
        rows = self._connection.execute(
            """
            SELECT amino_acid_symbol
            FROM amino_acid_class_members
            WHERE class_name = ?
            """,
            (class_name,),
        )
        return frozenset(str(row[0]) for row in rows)

    def codon_table(self, molecule_type: str) -> dict[str, str]:
        rows = self._connection.execute(
            """
            SELECT sequence, amino_acid_symbol
            FROM codons
            WHERE genetic_code_id = 1 AND molecule_type = ?
            """,
            (molecule_type,),
        )
        return {str(row[0]): str(row[1]) for row in rows}

    def stop_codons(self, molecule_type: str) -> frozenset[str]:
        rows = self._connection.execute(
            """
            SELECT sequence
            FROM codons
            WHERE genetic_code_id = 1
              AND molecule_type = ?
              AND is_stop = 1
            """,
            (molecule_type,),
        )
        return frozenset(str(row[0]) for row in rows)

    def start_codon(self, molecule_type: str) -> str:
        row = self._connection.execute(
            """
            SELECT sequence
            FROM codons
            WHERE genetic_code_id = 1
              AND molecule_type = ?
              AND is_start = 1
            """,
            (molecule_type,),
        ).fetchone()
        if row is None:
            raise RuntimeError(f"No start codon found for {molecule_type}.")
        return str(row[0])

    def constant(self, key: str) -> float:
        row = self._connection.execute(
            "SELECT value FROM reference_constants WHERE key = ?",
            (key,),
        ).fetchone()
        if row is None:
            raise RuntimeError(f"Missing biochemical reference constant: {key}.")
        return float(row[0])

    def close(self) -> None:
        self._connection.close()
        os.unlink(self._database_path)

    def __enter__(self) -> "ReferenceData":
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


@lru_cache(maxsize=1)
def get_reference_data() -> ReferenceData:
    return ReferenceData()
