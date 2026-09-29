"""Remove biochemical reference data from the API database.

Revision ID: 0004_remove_reference_data
Revises: 0003_merge_heads
"""

from typing import Sequence

from alembic import op


revision: str = "0004_remove_reference_data"
down_revision: str | Sequence[str] | None = "0003_merge_heads"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.drop_table("reference_constants")
    op.drop_table("codon_usage")
    op.drop_table("organisms")
    op.drop_table("codons")
    op.drop_table("nucleotides")
    op.drop_table("genetic_codes")
    op.drop_table("amino_acid_class_members")
    op.drop_table("amino_acid_classes")
    op.drop_table("reference_sources")
    op.drop_table("amino_acids")


def downgrade() -> None:
    raise RuntimeError(
        "Biochemical reference data is now owned by the BioSTAR engine as bundled Python data."
    )
