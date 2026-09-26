from __future__ import annotations

from decimal import Decimal

from sqlalchemy import Boolean, ForeignKey, Integer, Numeric, String, Text, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship

from .base import Base


class AminoAcid(Base):
    __tablename__ = "amino_acids"

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    symbol: Mapped[str] = mapped_column(String(1), unique=True, nullable=False)
    abbreviation: Mapped[str] = mapped_column(String(4), nullable=False)
    name: Mapped[str] = mapped_column(String(100), nullable=False)
    molecular_weight: Mapped[Decimal] = mapped_column(Numeric(10, 4), nullable=False)
    hydrophobicity: Mapped[Decimal] = mapped_column(Numeric(6, 3), nullable=False)
    alpha_helix: Mapped[Decimal] = mapped_column(Numeric(6, 3), nullable=False)
    beta_sheet: Mapped[Decimal] = mapped_column(Numeric(6, 3), nullable=False)
    pka: Mapped[Decimal | None] = mapped_column(Numeric(6, 3))
    pkb: Mapped[Decimal | None] = mapped_column(Numeric(6, 3))
    pkr: Mapped[Decimal | None] = mapped_column(Numeric(6, 3))

    classes: Mapped[list[AminoAcidClassMember]] = relationship(back_populates="amino_acid", cascade="all, delete-orphan")


class AminoAcidClass(Base):
    __tablename__ = "amino_acid_classes"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String(50), unique=True, nullable=False)
    description: Mapped[str | None] = mapped_column(Text)
    members: Mapped[list[AminoAcidClassMember]] = relationship(back_populates="amino_acid_class", cascade="all, delete-orphan")


class AminoAcidClassMember(Base):
    __tablename__ = "amino_acid_class_members"
    amino_acid_id: Mapped[int] = mapped_column(ForeignKey("amino_acids.id", ondelete="CASCADE"), primary_key=True)
    class_id: Mapped[int] = mapped_column(ForeignKey("amino_acid_classes.id", ondelete="CASCADE"), primary_key=True)
    amino_acid: Mapped[AminoAcid] = relationship(back_populates="classes")
    amino_acid_class: Mapped[AminoAcidClass] = relationship(back_populates="members")


class GeneticCode(Base):
    __tablename__ = "genetic_codes"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String(100), unique=True, nullable=False)
    ncbi_id: Mapped[int | None] = mapped_column(Integer, unique=True)
    description: Mapped[str | None] = mapped_column(Text)
    codons: Mapped[list[Codon]] = relationship(back_populates="genetic_code", cascade="all, delete-orphan")


class Nucleotide(Base):
    __tablename__ = "nucleotides"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    symbol: Mapped[str] = mapped_column(String(1), unique=True, nullable=False)
    name: Mapped[str] = mapped_column(String(50), nullable=False)
    molecule_type: Mapped[str] = mapped_column(String(10), nullable=False)


class Codon(Base):
    __tablename__ = "codons"
    __table_args__ = (UniqueConstraint("genetic_code_id", "molecule_type", "sequence", name="uq_codon_code_type_sequence"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    sequence: Mapped[str] = mapped_column(String(3), nullable=False)
    molecule_type: Mapped[str] = mapped_column(String(3), nullable=False)
    genetic_code_id: Mapped[int] = mapped_column(ForeignKey("genetic_codes.id", ondelete="CASCADE"), nullable=False)
    amino_acid_id: Mapped[int | None] = mapped_column(ForeignKey("amino_acids.id", ondelete="RESTRICT"))
    is_start: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    is_stop: Mapped[bool] = mapped_column(Boolean, nullable=False, default=False)
    genetic_code: Mapped[GeneticCode] = relationship(back_populates="codons")
    amino_acid: Mapped[AminoAcid | None] = relationship()


class Organism(Base):
    __tablename__ = "organisms"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    scientific_name: Mapped[str] = mapped_column(String(255), unique=True, nullable=False)
    common_name: Mapped[str | None] = mapped_column(String(255))
    taxonomy_id: Mapped[int | None] = mapped_column(Integer, unique=True)
    codon_usage: Mapped[list[CodonUsage]] = relationship(back_populates="organism", cascade="all, delete-orphan")


class ReferenceSource(Base):
    __tablename__ = "reference_sources"
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    name: Mapped[str] = mapped_column(String(255), unique=True, nullable=False)
    version: Mapped[str | None] = mapped_column(String(100))
    url: Mapped[str | None] = mapped_column(String(500))
    description: Mapped[str | None] = mapped_column(Text)


class CodonUsage(Base):
    __tablename__ = "codon_usage"
    __table_args__ = (UniqueConstraint("organism_id", "codon_id", "source_id", name="uq_codon_usage_organism_codon_source"),)
    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    organism_id: Mapped[int] = mapped_column(ForeignKey("organisms.id", ondelete="CASCADE"), nullable=False)
    codon_id: Mapped[int] = mapped_column(ForeignKey("codons.id", ondelete="CASCADE"), nullable=False)
    source_id: Mapped[int | None] = mapped_column(ForeignKey("reference_sources.id", ondelete="SET NULL"))
    frequency: Mapped[Decimal] = mapped_column(Numeric(12, 8), nullable=False)
    count: Mapped[int | None] = mapped_column(Integer)
    organism: Mapped[Organism] = relationship(back_populates="codon_usage")
    codon: Mapped[Codon] = relationship()
    source: Mapped[ReferenceSource | None] = relationship()
