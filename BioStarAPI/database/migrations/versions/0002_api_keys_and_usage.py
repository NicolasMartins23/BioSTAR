"""Add API key authentication and request quota tables.

Revision ID: 0002_api_keys_and_usage
Revises: 0001_initial_schema
"""

from typing import Sequence

from alembic import op
import sqlalchemy as sa


revision: str = "0002_api_keys_and_usage"
down_revision: str | Sequence[str] | None = "0001_initial_schema"
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    op.create_table(
        "api_keys",
        sa.Column("id", sa.Integer(), primary_key=True),
        sa.Column("name", sa.String(100)),
        sa.Column("key_hash", sa.String(64), nullable=False),
        sa.Column("key_prefix", sa.String(16), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False, server_default=sa.func.now()),
        sa.Column("revoked_at", sa.DateTime(timezone=True)),
        sa.UniqueConstraint("key_hash", name="uq_api_keys_key_hash"),
    )

    op.create_table(
        "api_daily_usage",
        sa.Column("usage_date", sa.Date(), nullable=False),
        sa.Column("scope", sa.String(20), nullable=False),
        sa.Column("scope_id", sa.String(100), nullable=False),
        sa.Column("request_count", sa.Integer(), nullable=False, server_default="0"),
        sa.PrimaryKeyConstraint("usage_date", "scope", "scope_id"),
    )


def downgrade() -> None:
    op.drop_table("api_daily_usage")
    op.drop_table("api_keys")
