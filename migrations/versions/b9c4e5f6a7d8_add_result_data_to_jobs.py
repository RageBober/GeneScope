"""Add result_data to jobs

Revision ID: b9c4e5f6a7d8
Revises: a8b65a5db1b4
Create Date: 2025-11-29 20:30:00.000000

"""
from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision: str = 'b9c4e5f6a7d8'
down_revision: Union[str, None] = 'a8b65a5db1b4'
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    """Add result_data JSON column to jobs table."""
    op.add_column(
        'jobs',
        sa.Column('result_data', sa.JSON(), nullable=True)
    )


def downgrade() -> None:
    """Remove result_data column from jobs table."""
    op.drop_column('jobs', 'result_data')
