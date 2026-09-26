from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from BioStar.data.biochemistry import (
    AMINOACID_TABLE,
    AMINOACIDS_AROMATIC,
    AMINOACIDS_NEGATIVE,
    AMINOACIDS_NONPOLAR,
    AMINOACIDS_POLAR,
    AMINOACIDS_POSITIVE,
    TABLE_DNA_CODON_TO_AMINOACID,
    TABLE_RNA_CODON_TO_AMINOACID,
)
from BioStarAPI.database.connection import engine
from BioStarAPI.database.models import (
    AminoAcid,
    AminoAcidClass,
    AminoAcidClassMember,
    Codon,
    GeneticCode,
    Nucleotide,
    ReferenceSource,
)


def seed_database() -> None:
    with Session(engine) as session:
        _seed_amino_acids(session)
        _seed_amino_acid_classes(session)
        _seed_nucleotides(session)
        _seed_genetic_code(session)
        session.commit()


def _seed_amino_acids(session: Session) -> None:
    existing: set[str] = set(
        session.scalars(select(AminoAcid.symbol))
    )

    for symbol, data in AMINOACID_TABLE.items():
        if symbol in existing:
            continue

        session.add(
            AminoAcid(
                symbol=symbol,
                abbreviation=data["abbreviation"],
                name=data["name"],
                molecular_weight=data["weight"],
                hydrophobicity=data["hydrophobicity"],
                alpha_helix=data["alpha_helix"],
                beta_sheet=data["beta_sheet"],
                pka=data.get("pKa"),
                pkb=data.get("pKb"),
                pkr=data.get("pKr"),
            )
        )


def _seed_amino_acid_classes(session: Session) -> None:
    class_members: dict[str, list[str]] = {
        "aromatic": AMINOACIDS_AROMATIC,
        "nonpolar": AMINOACIDS_NONPOLAR,
        "polar": AMINOACIDS_POLAR,
        "positive": AMINOACIDS_POSITIVE,
        "negative": AMINOACIDS_NEGATIVE,
    }

    amino_acids = {
        amino_acid.symbol: amino_acid
        for amino_acid in session.scalars(select(AminoAcid))
    }

    for class_name, symbols in class_members.items():
        amino_acid_class = session.scalar(
            select(AminoAcidClass).where(AminoAcidClass.name == class_name)
        )

        if amino_acid_class is None:
            amino_acid_class = AminoAcidClass(name=class_name)
            session.add(amino_acid_class)
            session.flush()

        for symbol in symbols:
            amino_acid = amino_acids[symbol]
            if not any(member.amino_acid_id == amino_acid.id for member in amino_acid_class.members):
                amino_acid_class.members.append(
                    AminoAcidClassMember(
                        amino_acid=amino_acid,
                        amino_acid_class=amino_acid_class,
                    )
                )


def _seed_nucleotides(session: Session) -> None:
    nucleotides = {
        ("A", "Adenine", "DNA"),
        ("C", "Cytosine", "DNA"),
        ("G", "Guanine", "DNA"),
        ("T", "Thymine", "DNA"),
        ("U", "Uracil", "RNA"),
    }

    existing: set[str] = set(session.scalars(select(Nucleotide.symbol)))

    for symbol, name, molecule_type in nucleotides:
        if symbol not in existing:
            session.add(
                Nucleotide(
                    symbol=symbol,
                    name=name,
                    molecule_type=molecule_type,
                )
            )


def _seed_genetic_code(session: Session) -> None:
    genetic_code = session.scalar(
        select(GeneticCode).where(GeneticCode.ncbi_id == 1)
    )

    if genetic_code is None:
        genetic_code = GeneticCode(
            name="Standard",
            ncbi_id=1,
            description="Standard genetic code.",
        )
        session.add(genetic_code)
        session.flush()

    amino_acids = {
        amino_acid.symbol: amino_acid
        for amino_acid in session.scalars(select(AminoAcid))
    }

    source = session.scalar(
        select(ReferenceSource).where(ReferenceSource.name == "BioSTAR legacy DataBiochemistry")
    )
    if source is None:
        source = ReferenceSource(
            name="BioSTAR legacy DataBiochemistry",
            version="0.2.0",
            description="Initial reference data migrated from BioSTAR's in-memory biochemical tables.",
        )
        session.add(source)

    for molecule_type, table in (
        ("DNA", TABLE_DNA_CODON_TO_AMINOACID),
        ("RNA", TABLE_RNA_CODON_TO_AMINOACID),
    ):
        for sequence, symbol in table.items():
            exists = session.scalar(
                select(Codon).where(
                    Codon.genetic_code_id == genetic_code.id,
                    Codon.molecule_type == molecule_type,
                    Codon.sequence == sequence,
                )
            )
            if exists is not None:
                continue

            session.add(
                Codon(
                    sequence=sequence,
                    molecule_type=molecule_type,
                    genetic_code=genetic_code,
                    amino_acid=amino_acids[symbol],
                    is_start=symbol == "M" and sequence in {"ATG", "AUG"},
                    is_stop=symbol == "*",
                )
            )


if __name__ == "__main__":
    seed_database()
