"""Create initial BioSTAR database schema.

Revision ID: 0001_initial_schema
Revises:
"""

from typing import Sequence

from alembic import op
import sqlalchemy as sa


revision: str = "0001_initial_schema"
down_revision: str | Sequence[str] | None = None
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.create_table(
        "amino_acids",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("symbol", sa.String(1), nullable=False),
        sa.Column("abbreviation", sa.String(3), nullable=False),
        sa.Column("name", sa.String(100), nullable=False),
        sa.Column("molecular_weight", sa.Numeric(10, 4), nullable=False),
        sa.Column("hydrophobicity", sa.Numeric(6, 3), nullable=False),
        sa.Column("alpha_helix", sa.Numeric(6, 3), nullable=False),
        sa.Column("beta_sheet", sa.Numeric(6, 3), nullable=False),
        sa.Column("pka", sa.Numeric(6, 3)),
        sa.Column("pkb", sa.Numeric(6, 3)),
        sa.Column("pkr", sa.Numeric(6, 3)),
        sa.UniqueConstraint("symbol", name="uq_amino_acids_symbol"),
    )

    op.create_table(
        "amino_acid_classes",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("name", sa.String(50), nullable=False),
        sa.Column("description", sa.Text()),
        sa.UniqueConstraint("name", name="uq_amino_acid_classes_name"),
    )

    op.create_table(
        "amino_acid_class_members",
        sa.Column("amino_acid_id", sa.Integer(), nullable=False),
        sa.Column("class_id", sa.Integer(), nullable=False),
        sa.ForeignKeyConstraint(["amino_acid_id"], ["amino_acids.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["class_id"], ["amino_acid_classes.id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("amino_acid_id", "class_id"),
    )

    op.create_table(
        "genetic_codes",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("name", sa.String(100), nullable=False),
        sa.Column("ncbi_id", sa.Integer()),
        sa.Column("description", sa.Text()),
        sa.UniqueConstraint("name", name="uq_genetic_codes_name"),
        sa.UniqueConstraint("ncbi_id", name="uq_genetic_codes_ncbi_id"),
    )

    op.create_table(
        "nucleotides",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("symbol", sa.String(1), nullable=False),
        sa.Column("name", sa.String(50), nullable=False),
        sa.Column("molecule_type", sa.String(10), nullable=False),
        sa.UniqueConstraint("symbol", name="uq_nucleotides_symbol"),
    )

    op.create_table(
        "codons",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("sequence", sa.String(3), nullable=False),
        sa.Column("molecule_type", sa.String(3), nullable=False),
        sa.Column("genetic_code_id", sa.Integer(), nullable=False),
        sa.Column("amino_acid_id", sa.Integer()),
        sa.Column("is_start", sa.Boolean(), nullable=False, server_default=sa.false()),
        sa.Column("is_stop", sa.Boolean(), nullable=False, server_default=sa.false()),
        sa.ForeignKeyConstraint(["genetic_code_id"], ["genetic_codes.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["amino_acid_id"], ["amino_acids.id"], ondelete="RESTRICT"),
        sa.UniqueConstraint(
            "genetic_code_id",
            "molecule_type",
            "sequence",
            name="uq_codon_code_type_sequence",
        ),
    )

    op.create_table(
        "organisms",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("scientific_name", sa.String(255), nullable=False),
        sa.Column("common_name", sa.String(255)),
        sa.Column("taxonomy_id", sa.Integer()),
        sa.UniqueConstraint("scientific_name", name="uq_organisms_scientific_name"),
        sa.UniqueConstraint("taxonomy_id", name="uq_organisms_taxonomy_id"),
    )

    op.create_table(
        "reference_sources",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("name", sa.String(255), nullable=False),
        sa.Column("version", sa.String(100)),
        sa.Column("url", sa.String(500)),
        sa.Column("description", sa.Text()),
        sa.UniqueConstraint("name", name="uq_reference_sources_name"),
    )

    op.create_table(
        "codon_usage",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("organism_id", sa.Integer(), nullable=False),
        sa.Column("codon_id", sa.Integer(), nullable=False),
        sa.Column("source_id", sa.Integer()),
        sa.Column("frequency", sa.Numeric(12, 8), nullable=False),
        sa.Column("count", sa.Integer()),
        sa.ForeignKeyConstraint(["organism_id"], ["organisms.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["codon_id"], ["codons.id"], ondelete="CASCADE"),
        sa.ForeignKeyConstraint(["source_id"], ["reference_sources.id"], ondelete="SET NULL"),
        sa.UniqueConstraint(
            "organism_id",
            "codon_id",
            "source_id",
            name="uq_codon_usage_organism_codon_source",
        ),
    )

    op.create_table(
        "users",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("email", sa.String(320), nullable=False),
        sa.Column("password_hash", sa.String(255), nullable=False),
        sa.Column("is_active", sa.Boolean(), nullable=False, server_default=sa.true()),
        sa.UniqueConstraint("email", name="uq_users_email"),
    )

    op.create_table(
        "refresh_tokens",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("user_id", sa.Integer(), nullable=False),
        sa.Column("token_hash", sa.String(255), nullable=False),
        sa.Column("expires_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("revoked_at", sa.DateTime(timezone=True)),
        sa.ForeignKeyConstraint(["user_id"], ["users.id"], ondelete="CASCADE"),
        sa.UniqueConstraint("token_hash", name="uq_refresh_tokens_token_hash"),
    )


def downgrade() -> None:
    op.drop_table("refresh_tokens")
    op.drop_table("users")
    op.drop_table("codon_usage")
    op.drop_table("reference_sources")
    op.drop_table("organisms")
    op.drop_table("codons")
    op.drop_table("nucleotides")
    op.drop_table("genetic_codes")
    op.drop_table("amino_acid_class_members")
    op.drop_table("amino_acid_classes")
    op.drop_table("amino_acids")
