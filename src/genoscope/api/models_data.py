from __future__ import annotations

from datetime import datetime
from typing import Any

from sqlmodel import JSON
from sqlmodel import Column
from sqlmodel import Field
from sqlmodel import SQLModel


class Dataset(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    original_name: str
    path_parquet: str
    n_rows: int
    n_cols: int
    columns: list[str] = Field(sa_column=Column(JSON))
    mapping: dict[str, str] = Field(default={}, sa_column=Column(JSON))
    created_at: datetime = Field(default_factory=datetime.utcnow)


class FilterRun(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    dataset_id: int = Field(index=True)
    params: dict[str, Any] = Field(sa_column=Column(JSON))
    preset: str | None = Field(default=None, index=True)
    total: int = 0
    created_at: datetime = Field(default_factory=datetime.utcnow)


class Report(SQLModel, table=True):
    id: int | None = Field(default=None, primary_key=True)
    filter_run_id: int = Field(index=True)
    path: str
    fmt: str = "html"
    created_at: datetime = Field(default_factory=datetime.utcnow)
