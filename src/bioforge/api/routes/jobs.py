"""Jobs routes."""

import logging
from fastapi import APIRouter, HTTPException, status
from pydantic import BaseModel
from sqlalchemy import select

from bioforge.api.deps import SessionDep, CurrentActiveUser, ArqPoolDep
from bioforge.models import Job, Sample, Project, JobStatus, JobType

router = APIRouter(prefix="/jobs", tags=["jobs"])
logger = logging.getLogger("bioforge.api.jobs")


# ─────────────────────────────────────────────────────────────────────────────
# Schemas
# ─────────────────────────────────────────────────────────────────────────────

class JobCreate(BaseModel):
    sample_id: int
    job_type: str = JobType.FULL_PIPELINE.value
    parameters: dict | None = None


class JobResponse(BaseModel):
    id: int
    job_type: str
    status: str
    progress: float
    current_step: str | None
    result_path: str | None
    error_message: str | None
    task_id: str | None
    sample_id: int

    model_config = {"from_attributes": True}


# Task function mapping
TASK_FUNCTIONS = {
    JobType.ALIGNMENT.value: "run_alignment",
    JobType.VARIANT_CALLING.value: "run_variant_calling",
    JobType.FULL_PIPELINE.value: "run_full_pipeline",
}


# ─────────────────────────────────────────────────────────────────────────────
# Routes
# ─────────────────────────────────────────────────────────────────────────────

@router.get("", response_model=list[JobResponse])
async def list_jobs(
    session: SessionDep,
    user: CurrentActiveUser,
    sample_id: int | None = None,
    status_filter: str | None = None,
):
    """List jobs for current user's samples."""
    query = (
        select(Job)
        .join(Sample)
        .join(Project)
        .where(Project.owner_id == user.id)
    )
    
    if sample_id:
        query = query.where(Job.sample_id == sample_id)
    if status_filter:
        query = query.where(Job.status == status_filter)
    
    query = query.order_by(Job.created_at.desc())
    result = await session.execute(query)
    return result.scalars().all()


@router.post("", response_model=JobResponse, status_code=status.HTTP_201_CREATED)
async def create_job(
    data: JobCreate,
    session: SessionDep,
    user: CurrentActiveUser,
    arq_pool: ArqPoolDep,
):
    """Create and enqueue a new pipeline job."""
    # Verify sample ownership
    result = await session.execute(
        select(Sample)
        .join(Project)
        .where(Sample.id == data.sample_id, Project.owner_id == user.id)
    )
    sample = result.scalar_one_or_none()
    if not sample:
        raise HTTPException(status_code=404, detail="Sample not found")
    
    # Validate job type
    valid_types = [t.value for t in JobType]
    if data.job_type not in valid_types:
        raise HTTPException(
            status_code=400,
            detail=f"Invalid job type. Must be one of: {valid_types}",
        )
    
    # Create job
    job = Job(
        sample_id=data.sample_id,
        job_type=data.job_type,
        parameters=data.parameters,
        status=JobStatus.PENDING.value,
    )
    session.add(job)
    await session.commit()
    await session.refresh(job)
    
    # Enqueue task to ARQ
    task_func = TASK_FUNCTIONS.get(data.job_type)
    if task_func:
        try:
            arq_job = await arq_pool.enqueue_job(task_func, job.id)
            job.task_id = arq_job.job_id
            job.status = JobStatus.QUEUED.value
            await session.commit()
            await session.refresh(job)
            logger.info(f"Job {job.id} enqueued as task {arq_job.job_id}")
        except Exception as e:
            logger.error(f"Failed to enqueue job {job.id}: {e}")
            # Job stays in PENDING state, can be retried
    
    return job


@router.get("/{job_id}", response_model=JobResponse)
async def get_job(
    job_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
):
    """Get job by ID."""
    result = await session.execute(
        select(Job)
        .join(Sample)
        .join(Project)
        .where(Job.id == job_id, Project.owner_id == user.id)
    )
    job = result.scalar_one_or_none()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    return job


@router.post("/{job_id}/retry", response_model=JobResponse)
async def retry_job(
    job_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
    arq_pool: ArqPoolDep,
):
    """Retry a failed job."""
    result = await session.execute(
        select(Job)
        .join(Sample)
        .join(Project)
        .where(Job.id == job_id, Project.owner_id == user.id)
    )
    job = result.scalar_one_or_none()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if job.status not in [JobStatus.FAILED.value, JobStatus.CANCELLED.value]:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot retry job with status: {job.status}",
        )
    
    # Reset job
    job.status = JobStatus.PENDING.value
    job.progress = 0
    job.current_step = None
    job.error_message = None
    job.started_at = None
    job.completed_at = None
    await session.commit()
    
    # Re-enqueue
    task_func = TASK_FUNCTIONS.get(job.job_type)
    if task_func:
        try:
            arq_job = await arq_pool.enqueue_job(task_func, job.id)
            job.task_id = arq_job.job_id
            job.status = JobStatus.QUEUED.value
            await session.commit()
        except Exception as e:
            logger.error(f"Failed to re-enqueue job {job.id}: {e}")
    
    await session.refresh(job)
    return job


@router.post("/{job_id}/cancel", response_model=JobResponse)
async def cancel_job(
    job_id: int,
    session: SessionDep,
    user: CurrentActiveUser,
    arq_pool: ArqPoolDep,
):
    """Cancel a running or pending job."""
    result = await session.execute(
        select(Job)
        .join(Sample)
        .join(Project)
        .where(Job.id == job_id, Project.owner_id == user.id)
    )
    job = result.scalar_one_or_none()
    if not job:
        raise HTTPException(status_code=404, detail="Job not found")
    
    if job.status not in [JobStatus.PENDING.value, JobStatus.QUEUED.value, JobStatus.RUNNING.value]:
        raise HTTPException(
            status_code=400,
            detail=f"Cannot cancel job with status: {job.status}",
        )
    
    # Try to abort ARQ task
    if job.task_id:
        try:
            await arq_pool.abort_job(job.task_id)
        except Exception as e:
            logger.warning(f"Could not abort ARQ task {job.task_id}: {e}")
    
    job.status = JobStatus.CANCELLED.value
    await session.commit()
    await session.refresh(job)
    return job
