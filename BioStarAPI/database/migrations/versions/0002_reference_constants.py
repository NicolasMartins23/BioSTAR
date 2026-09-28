"""Add biochemical reference constants.

Revision ID: 0002_reference_constants
Revises: 0001_initial_schema
"""

from typing import Sequence

from alembic import op
import sqlalchemy as sa


revision: str = "0002_reference_constants"
down_revision: str | Sequence[str] | None = "0001_initial_schema"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.create_table(
        "reference_constants",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("key", sa.String(100), nullable=False),
        sa.Column("value", sa.Numeric(12, 6), nullable=False),
        sa.Column("description", sa.Text()),
        sa.UniqueConstraint("key", name="uq_reference_constants_key"),
    )


def downgrade() -> None:
    op.drop_table("reference_constants")
