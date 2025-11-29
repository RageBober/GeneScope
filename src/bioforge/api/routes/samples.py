"""Samples routes."""

from fastapi import APIRouter, HTTPException, UploadFile, File, status
from pydantic import BaseModel
from sqlalchemy import select

from bioforge.api.deps import SessionDep, CurrentActiveUser
from bioforge.config import settings
from bioforge.core.security import secure_filename, async_save_upload, validate_file
from bioforge.core.exceptions import FileTooLargeError, InvalidFileTypeError
from bioforge.models import Sample, Project, SampleStatus

router = APIRouter(prefix="/samples", tags=["samples"])


# ─────────────────────────────────────────────────────────────────────────────
# Schemas
# ─────────────────────────────────────────────────────────────────────────────

class SampleCreate(BaseModel):
    name: str
    project_id: int
    description: str | None = None
    sample_type: str | None = None
    organism: str = "Homo sapiens"
    tissue: str | None = None


class SampleResponse(BaseModel):
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
    project_id: int

    model_config = {"from_attributes": True}


# ─────────────────────────────────────────────────────────────────────────────
# Routes
# ─────────────────────────────────────────────────────────────────────────────

@router.get("", response_model=list[SampleResponse])
async def list_samples(
    session: SessionDep,
    user: CurrentActiveUser,
    project_id: int | None = None,
):
    """List samples (optionally filtered by project)."""
    query = (
        select(Sample)
        .join(Project)
        .where(Project.owner_id == user.id)
    )
    if project_id:
        query = query.where(Sample.project_id == project_id)
    
    result = await session.execute(query)
    return result.scalars().all()


@router.post("", response_model=SampleResponse, status_code=status.HTTP_201_CREATED)
async def create_sample(
    data: SampleCreate,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Create a new sample."""
    # Verify project ownership
    result = await session.execute(
        select(Project).where(
            Project.id == data.project_id,
            Project.owner_id == user.id,
        )
    )
    if not result.scalar_one_or_none():
        raise HTTPException(status_code=404, detail="Project not found")
    
    sample = Sample(
        name=data.name,
        description=data.description,
        sample_type=data.sample_type,
        organism=data.organism,
        tissue=data.tissue,
        project_id=data.project_id,
    )
    session.add(sample)
    await session.commit()
    await session.refresh(sample)
    return sample


@router.get("/{sample_id}", response_model=SampleResponse)
async def get_sample(
    sample_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Get sample by ID."""
    result = await session.execute(
        select(Sample)
        .join(Project)
        .where(Sample.id == sample_id, Project.owner_id == user.id)
    )
    sample = result.scalar_one_or_none()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    return sample


@router.post("/{sample_id}/upload")
async def upload_fastq(
    sample_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
    r1: UploadFile = File(...),
    r2: UploadFile | None = File(None),
):
    """Upload FASTQ files for a sample."""
    # Get sample
    result = await session.execute(
        select(Sample)
        .join(Project)
        .where(Sample.id == sample_id, Project.owner_id == user.id)
    )
    sample = result.scalar_one_or_none()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    
    # Save files
    sample_dir = settings.storage.samples_dir / str(sample.project_id) / str(sample.id)
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Save R1
    r1_filename = secure_filename(r1.filename or "r1.fastq.gz")
    r1_path = sample_dir / r1_filename
    try:
        await async_save_upload(r1, r1_path, settings.storage.max_file_size)
        validate_file(r1_path, settings.storage.allowed_extensions, settings.storage.max_file_size)
    except (FileTooLargeError, InvalidFileTypeError) as e:
        r1_path.unlink(missing_ok=True)
        raise HTTPException(status_code=e.status_code, detail=e.message)
    
    sample.fastq_r1 = str(r1_path)
    
    # Save R2 if provided
    if r2:
        r2_filename = secure_filename(r2.filename or "r2.fastq.gz")
        r2_path = sample_dir / r2_filename
        try:
            await async_save_upload(r2, r2_path, settings.storage.max_file_size)
            validate_file(r2_path, settings.storage.allowed_extensions, settings.storage.max_file_size)
        except (FileTooLargeError, InvalidFileTypeError) as e:
            r2_path.unlink(missing_ok=True)
            raise HTTPException(status_code=e.status_code, detail=e.message)
        sample.fastq_r2 = str(r2_path)
    
    sample.status = SampleStatus.UPLOADED.value
    await session.commit()
    
    return {"status": "uploaded", "sample_id": sample.id}


@router.delete("/{sample_id}", status_code=status.HTTP_204_NO_CONTENT)
async def delete_sample(
    sample_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Delete sample."""
    result = await session.execute(
        select(Sample)
        .join(Project)
        .where(Sample.id == sample_id, Project.owner_id == user.id)
    )
    sample = result.scalar_one_or_none()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    
    await session.delete(sample)
    await session.commit()
