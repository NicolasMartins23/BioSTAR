"""Expand amino acid abbreviations to four characters.

Revision ID: 0002_expand_amino_acid_abbreviation
Revises: 0001_initial_schema
"""

from typing import Sequence

from alembic import op


revision: str = "0002_expand_amino_acid_abbreviation"
down_revision: str | Sequence[str] | None = "0001_initial_schema"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.alter_column(
        "amino_acids",
        "abbreviation",
        existing_type=__import__("sqlalchemy").String(3),
        type_=__import__("sqlalchemy").String(4),
        existing_nullable=False,
    )


def downgrade() -> None:
    op.alter_column(
        "amino_acids",
        "abbreviation",
        existing_type=__import__("sqlalchemy").String(4),
        type_=__import__("sqlalchemy").String(3),
        existing_nullable=False,
    )
