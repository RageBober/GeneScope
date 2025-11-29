"""Variants routes."""

import logging
from typing import Literal

from fastapi import APIRouter, HTTPException, Query
from pydantic import BaseModel, Field
from sqlalchemy import select, func, and_

from bioforge.api.deps import SessionDep, CurrentActiveUser
from bioforge.models import Variant, Sample, Project

router = APIRouter(prefix="/variants", tags=["variants"])
logger = logging.getLogger("bioforge.api.variants")


# ─────────────────────────────────────────────────────────────────────────────
# Schemas
# ─────────────────────────────────────────────────────────────────────────────

class VariantResponse(BaseModel):
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
    rare: int


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
# Routes
# ─────────────────────────────────────────────────────────────────────────────

@router.get("/sample/{sample_id}", response_model=list[VariantResponse])
async def list_variants(
    sample_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
    chromosome: str | None = None,
    gene: str | None = None,
    impact: str | None = Query(None, description="Comma-separated: HIGH,MODERATE,LOW"),
    min_quality: float | None = None,
    max_gnomad_af: float | None = None,
    clinvar_only: bool = False,
    pathogenic_only: bool = False,
    limit: int = Query(100, le=1000),
    offset: int = 0,
):
    """List variants for a sample with filtering."""
    # Verify sample ownership
    result = await session.execute(
        select(Sample)
        .join(Project)
        .where(Sample.id == sample_id, Project.owner_id == user.id)
    )
    sample = result.scalar_one_or_none()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    
    # Build query
    query = select(Variant).where(Variant.sample_id == sample_id)
    
    if chromosome:
        query = query.where(Variant.chromosome == chromosome)
    
    if gene:
        query = query.where(Variant.gene_symbol.ilike(f"%{gene}%"))
    
    if impact:
        impacts = [i.strip().upper() for i in impact.split(",")]
        query = query.where(Variant.impact.in_(impacts))
    
    if min_quality is not None:
        query = query.where(Variant.quality >= min_quality)
    
    if max_gnomad_af is not None:
        query = query.where(
            (Variant.gnomad_af == None) | (Variant.gnomad_af <= max_gnomad_af)
        )
    
    if clinvar_only:
        query = query.where(Variant.clinvar_id != None)
    
    if pathogenic_only:
        query = query.where(
            Variant.clinvar_significance.ilike("%pathogenic%")
        )
    
    # Order and paginate
    query = query.order_by(Variant.chromosome, Variant.position)
    query = query.offset(offset).limit(limit)
    
    result = await session.execute(query)
    return result.scalars().all()


@router.get("/sample/{sample_id}/stats", response_model=VariantStats)
async def get_variant_stats(
    sample_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Get variant statistics for a sample."""
    # Verify sample ownership
    result = await session.execute(
        select(Sample)
        .join(Project)
        .where(Sample.id == sample_id, Project.owner_id == user.id)
    )
    if not result.scalar_one_or_none():
        raise HTTPException(status_code=404, detail="Sample not found")
    
    # Count total
    total = await session.scalar(
        select(func.count()).where(Variant.sample_id == sample_id)
    )
    
    # Count by type
    snv = await session.scalar(
        select(func.count()).where(
            Variant.sample_id == sample_id,
            Variant.variant_type == "SNV"
        )
    )
    
    indel = await session.scalar(
        select(func.count()).where(
            Variant.sample_id == sample_id,
            Variant.variant_type == "INDEL"
        )
    )
    
    # High impact
    high_impact = await session.scalar(
        select(func.count()).where(
            Variant.sample_id == sample_id,
            Variant.impact == "HIGH"
        )
    )
    
    # Pathogenic
    pathogenic = await session.scalar(
        select(func.count()).where(
            Variant.sample_id == sample_id,
            Variant.clinvar_significance.ilike("%pathogenic%")
        )
    )
    
    # Rare (AF < 1%)
    rare = await session.scalar(
        select(func.count()).where(
            Variant.sample_id == sample_id,
            (Variant.gnomad_af == None) | (Variant.gnomad_af < 0.01)
        )
    )
    
    return VariantStats(
        total=total or 0,
        snv=snv or 0,
        indel=indel or 0,
        high_impact=high_impact or 0,
        pathogenic=pathogenic or 0,
        rare=rare or 0,
    )


@router.get("/{variant_id}", response_model=VariantDetail)
async def get_variant(
    variant_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Get variant by ID with full details."""
    result = await session.execute(
        select(Variant)
        .join(Sample)
        .join(Project)
        .where(Variant.id == variant_id, Project.owner_id == user.id)
    )
    variant = result.scalar_one_or_none()
    if not variant:
        raise HTTPException(status_code=404, detail="Variant not found")
    return variant


@router.get("/search", response_model=list[VariantResponse])
async def search_variants(
    session: SessionDep,
    user: CurrentActiveUser,
    query: str = Query(..., min_length=2, description="Search by gene, rsID, or position"),
    project_id: int | None = None,
    limit: int = Query(50, le=200),
):
    """Search variants across user's samples."""
    # Base query with ownership filter
    stmt = (
        select(Variant)
        .join(Sample)
        .join(Project)
        .where(Project.owner_id == user.id)
    )
    
    if project_id:
        stmt = stmt.where(Project.id == project_id)
    
    # Search logic
    if query.startswith("rs"):
        # rsID search
        stmt = stmt.where(Variant.clinvar_id == query)
    elif ":" in query:
        # Position search (chr:pos)
        parts = query.split(":")
        if len(parts) == 2:
            chrom = parts[0].replace("chr", "")
            try:
                pos = int(parts[1])
                stmt = stmt.where(
                    Variant.chromosome == chrom,
                    Variant.position == pos
                )
            except ValueError:
                pass
    else:
        # Gene search
        stmt = stmt.where(Variant.gene_symbol.ilike(f"%{query}%"))
    
    stmt = stmt.limit(limit)
    
    result = await session.execute(stmt)
    return result.scalars().all()


@router.get("/genes/top", response_model=list[dict])
async def top_genes_by_variants(
    session: SessionDep,
    user: CurrentActiveUser,
    sample_id: int | None = None,
    project_id: int | None = None,
    impact: str | None = Query(None, description="Filter by impact: HIGH,MODERATE"),
    limit: int = Query(20, le=100),
):
    """Get top genes by variant count."""
    stmt = (
        select(
            Variant.gene_symbol,
            func.count(Variant.id).label("variant_count")
        )
        .join(Sample)
        .join(Project)
        .where(
            Project.owner_id == user.id,
            Variant.gene_symbol != None
        )
    )
    
    if sample_id:
        stmt = stmt.where(Variant.sample_id == sample_id)
    
    if project_id:
        stmt = stmt.where(Project.id == project_id)
    
    if impact:
        impacts = [i.strip().upper() for i in impact.split(",")]
        stmt = stmt.where(Variant.impact.in_(impacts))
    
    stmt = (
        stmt
        .group_by(Variant.gene_symbol)
        .order_by(func.count(Variant.id).desc())
        .limit(limit)
    )
    
    result = await session.execute(stmt)
    return [
        {"gene": row[0], "count": row[1]}
        for row in result.all()
    ]
