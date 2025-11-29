"""Pydantic schemas for API request/response models.

Re-exports schemas from route modules for convenience.
"""

from pydantic import BaseModel, EmailStr, Field
from datetime import datetime
from typing import Any


# ─────────────────────────────────────────────────────────────────────────────
# Auth Schemas
# ─────────────────────────────────────────────────────────────────────────────

class UserCreate(BaseModel):
    """Schema for user registration."""
    email: EmailStr
    password: str = Field(..., min_length=8)
    full_name: str | None = None


class UserLogin(BaseModel):
    """Schema for user login."""
    email: EmailStr
    password: str


class Token(BaseModel):
    """JWT token response."""
    access_token: str
    refresh_token: str
    token_type: str = "bearer"


class TokenPayload(BaseModel):
    """JWT token payload."""
    sub: str
    exp: datetime
    type: str


class UserResponse(BaseModel):
    """User response schema."""
    id: int
    email: str
    full_name: str | None
    is_active: bool
    is_verified: bool
    created_at: datetime

    model_config = {"from_attributes": True}


# ─────────────────────────────────────────────────────────────────────────────
# Project Schemas
# ─────────────────────────────────────────────────────────────────────────────

class ProjectCreate(BaseModel):
    """Schema for creating a project."""
    name: str = Field(..., min_length=1, max_length=255)
    description: str | None = None
    reference_genome: str = "GRCh38"


class ProjectUpdate(BaseModel):
    """Schema for updating a project."""
    name: str | None = Field(None, min_length=1, max_length=255)
    description: str | None = None
    reference_genome: str | None = None


class ProjectResponse(BaseModel):
    """Project response schema."""
    id: int
    name: str
    description: str | None
    reference_genome: str
    owner_id: int
    created_at: datetime
    updated_at: datetime

    model_config = {"from_attributes": True}


# ─────────────────────────────────────────────────────────────────────────────
# Sample Schemas
# ─────────────────────────────────────────────────────────────────────────────

class SampleCreate(BaseModel):
    """Schema for creating a sample."""
    name: str = Field(..., min_length=1, max_length=255)
    project_id: int
    description: str | None = None
    sample_type: str | None = None
    organism: str = "Homo sapiens"
    tissue: str | None = None


class SampleUpdate(BaseModel):
    """Schema for updating a sample."""
    name: str | None = Field(None, min_length=1, max_length=255)
    description: str | None = None
    sample_type: str | None = None
    tissue: str | None = None


class SampleResponse(BaseModel):
    """Sample response schema."""
    id: int
    name: str
    description: str | None
    sample_type: str | None
    organism: str
    tissue: str | None
    status: str
    fastq_r1: str | None
    fastq_r2: str | None
    bam_path: str | None
    vcf_path: str | None
    qc_metrics: dict | None
    project_id: int
    created_at: datetime
    updated_at: datetime

    model_config = {"from_attributes": True}


# ─────────────────────────────────────────────────────────────────────────────
# Job Schemas
# ─────────────────────────────────────────────────────────────────────────────

class JobCreate(BaseModel):
    """Schema for creating a job."""
    sample_id: int
    job_type: str = "full_pipeline"
    parameters: dict | None = None


class JobResponse(BaseModel):
    """Job response schema."""
    id: int
    job_type: str
    status: str
    progress: float
    current_step: str | None
    result_path: str | None
    error_message: str | None
    task_id: str | None
    parameters: dict | None
    sample_id: int
    created_at: datetime
    started_at: datetime | None
    completed_at: datetime | None

    model_config = {"from_attributes": True}


# ─────────────────────────────────────────────────────────────────────────────
# Variant Schemas
# ─────────────────────────────────────────────────────────────────────────────

class VariantResponse(BaseModel):
    """Variant response schema."""
    id: int
    chromosome: str
    position: int
    ref: str
    alt: str
    variant_type: str
    quality: float | None
    filter_status: str
    genotype: str | None
    depth: int | None
    gene_symbol: str | None
    consequence: str | None
    impact: str | None
    hgvs_c: str | None
    hgvs_p: str | None
    gnomad_af: float | None
    clinvar_significance: str | None
    cadd_score: float | None
    sample_id: int

    model_config = {"from_attributes": True}


class VariantDetail(VariantResponse):
    """Full variant details."""
    gene_id: str | None
    transcript_id: str | None
    gnomad_af_popmax: float | None
    clinvar_id: str | None
    clinvar_conditions: str | None
    revel_score: float | None
    sift_pred: str | None
    polyphen_pred: str | None
    allele_depth: str | None
    annotations: dict | None


class VariantStats(BaseModel):
    """Variant statistics for a sample."""
    total: int
    snv: int
    indel: int
    high_impact: int
    pathogenic: int
    rare: int  # AF < 1%


class VariantFilter(BaseModel):
    """Filter parameters for variant queries."""
    chromosome: str | None = None
    gene: str | None = None
    impact: list[str] | None = None
    min_quality: float | None = None
    max_gnomad_af: float | None = None
    clinvar_only: bool = False
    pathogenic_only: bool = False


# ─────────────────────────────────────────────────────────────────────────────
# Common Schemas
# ─────────────────────────────────────────────────────────────────────────────

class PaginatedResponse(BaseModel):
    """Generic paginated response."""
    items: list[Any]
    total: int
    page: int
    page_size: int
    pages: int


class HealthResponse(BaseModel):
    """Health check response."""
    status: str
    version: str
    environment: str


class MessageResponse(BaseModel):
    """Simple message response."""
    message: str
    detail: str | None = None


# ─────────────────────────────────────────────────────────────────────────────
# Exports
# ─────────────────────────────────────────────────────────────────────────────

__all__ = [
    # Auth
    "UserCreate",
    "UserLogin",
    "Token",
    "TokenPayload",
    "UserResponse",
    # Project
    "ProjectCreate",
    "ProjectUpdate",
    "ProjectResponse",
    # Sample
    "SampleCreate",
    "SampleUpdate",
    "SampleResponse",
    # Job
    "JobCreate",
    "JobResponse",
    # Variant
    "VariantResponse",
    "VariantDetail",
    "VariantStats",
    "VariantFilter",
    # Common
    "PaginatedResponse",
    "HealthResponse",
    "MessageResponse",
]
