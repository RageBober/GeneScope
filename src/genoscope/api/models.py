# src/genoscope/api/models.py
from datetime import datetime

from sqlmodel import Field
from sqlmodel import SQLModel


class Job(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    filename: str
    status: str = "PENDING"  # PENDING / RUNNING / DONE / ERROR
    created_at: datetime = Field(default_factory=datetime.utcnow)
    updated_at: datetime = Field(default_factory=datetime.utcnow)
    result_path: str | None = None
