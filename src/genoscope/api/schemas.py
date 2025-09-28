from __future__ import annotations

from typing import Any

from pydantic import BaseModel


class UploadResponse(BaseModel):
    dataset_id: int
    rows: int
    cols: int
    columns: list[str]
    head: list[dict[str, Any]]
    mapping: dict[str, str]


class EvidenceCard(BaseModel):
    name: str  # правило/пресет, например "Rare monogenic"
    source: str  # gnomAD, ClinVar и т.д.
    version: str  # "v4.0"
    date: str  # ISO 8601
    details: dict[str, Any] = {}


class FilterParams(BaseModel):
    dataset_id: int
    preset: str | None = None
    min_af: float | None = None
    max_af: float | None = None
    min_qual: float | None = None
    genes: list[str] | None = None
    chroms: list[str] | None = None
    limit: int = 50
    create_report: bool = False
    locale: str = "ru"  # ru/kz/en
    clinic: str | None = None  # название клиники для отчёта
    logo_url: str | None = None


class FilterResponse(BaseModel):
    filter_run_id: int
    total: int
    preview: list[dict[str, Any]]
    applied: list[EvidenceCard]
    report_id: int | None = None
