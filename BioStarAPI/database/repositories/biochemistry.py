from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from BioStar.engine.biochemistry import AminoAcidData, BiochemistryData
from BioStarAPI.database.models import AminoAcid, AminoAcidClass, Codon, GeneticCode
from BioStarAPI.resources.exceptions import BioStarAPIError
from BioStarAPI.resources.messages import MessageCode


class BiochemistryRepository:
    """Loads immutable biochemical reference data from PostgreSQL."""

    def __init__(self, session: Session) -> None:
        self.session = session

    def get_standard_data(self) -> BiochemistryData:
        genetic_code = self.session.scalar(
            select(GeneticCode).where(GeneticCode.ncbi_id == 1)
        )
        if genetic_code is None:
            raise BioStarAPIError(500, MessageCode.STANDARD_GENETIC_CODE_NOT_SEEDED)

        amino_acids = list(self.session.scalars(select(AminoAcid)))
        if not amino_acids:
            raise BioStarAPIError(500, MessageCode.AMINO_ACID_DATA_NOT_SEEDED)

        classes = list(
            self.session.scalars(
                select(AminoAcidClass).where(
                    AminoAcidClass.name.in_(
                        ["aromatic", "nonpolar", "polar", "positive", "negative"]
                    )
                )
            )
        )
        class_symbols = {
            item.name: {
                member.amino_acid.symbol
                for member in item.members
                if member.amino_acid is not None
            }
            for item in classes
        }

        codons = list(
            self.session.scalars(
                select(Codon).where(Codon.genetic_code_id == genetic_code.id)
            )
        )
        dna_codons = {
            codon.sequence: codon.amino_acid.symbol
            for codon in codons
            if codon.molecule_type == "DNA" and codon.amino_acid is not None
        }
        rna_codons = {
            codon.sequence: codon.amino_acid.symbol
            for codon in codons
            if codon.molecule_type == "RNA" and codon.amino_acid is not None
        }

        return BiochemistryData(
            amino_acids={
                amino_acid.symbol: AminoAcidData(
                    symbol=amino_acid.symbol,
                    molecular_weight=float(amino_acid.molecular_weight),
                    hydrophobicity=float(amino_acid.hydrophobicity),
                    alpha_helix=float(amino_acid.alpha_helix),
                    beta_sheet=float(amino_acid.beta_sheet),
                    pkr=float(amino_acid.pkr) if amino_acid.pkr is not None else None,
                )
                for amino_acid in amino_acids
                if amino_acid.symbol != "*"
            },
            aromatic=frozenset(class_symbols.get("aromatic", set())),
            nonpolar=frozenset(class_symbols.get("nonpolar", set())),
            polar=frozenset(class_symbols.get("polar", set())),
            positive=frozenset(class_symbols.get("positive", set())),
            negative=frozenset(class_symbols.get("negative", set())),
            dna_codons=dna_codons,
            rna_codons=rna_codons,
            dna_stop_codons=frozenset(
                codon.sequence for codon in codons
                if codon.molecule_type == "DNA" and codon.is_stop
            ),
            rna_stop_codons=frozenset(
                codon.sequence for codon in codons
                if codon.molecule_type == "RNA" and codon.is_stop
            ),
            dna_start_codon=next(
                codon.sequence for codon in codons
                if codon.molecule_type == "DNA" and codon.is_start
            ),
            rna_start_codon=next(
                codon.sequence for codon in codons
                if codon.molecule_type == "RNA" and codon.is_start
            ),
        )
