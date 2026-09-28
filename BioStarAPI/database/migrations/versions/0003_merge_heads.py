"""Merge API key and reference constant migration branches.

Revision ID: 0003_merge_heads
Revises: 0002_api_keys_and_usage, 0002_reference_constants
"""

from typing import Sequence

revision: str = "0003_merge_heads"
down_revision: str | Sequence[str] | None = (
    "0002_api_keys_and_usage",
    "0002_reference_constants",
)
branch_labels: str | Sequence[str] | None = None
depends_on: str | Sequence[str] | None = None


def upgrade() -> None:
    pass


def downgrade() -> None:
    pass
