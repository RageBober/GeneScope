"""Pipeline orchestration service."""

import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from sqlalchemy import select
from sqlalchemy.ext.asyncio import AsyncSession

from bioforge.config import settings
from bioforge.models import Job, Sample, JobStatus, JobType
from bioforge.services.storage import storage_service

logger = logging.getLogger("bioforge.pipeline")


class PipelineService:
    """Service for orchestrating genomic analysis pipelines."""
    
    def __init__(self, session: AsyncSession):
        self.session = session
    
    async def create_job(
        self,
        sample_id: int,
        job_type: str,
        parameters: dict[str, Any] | None = None,
    ) -> Job:
        """Create a new pipeline job."""
        job = Job(
            sample_id=sample_id,
            job_type=job_type,
            parameters=parameters or {},
            status=JobStatus.PENDING.value,
        )
        self.session.add(job)
        await self.session.commit()
        await self.session.refresh(job)
        
        logger.info(f"Created job {job.id} for sample {sample_id}, type: {job_type}")
        return job
    
    async def update_job_status(
        self,
        job_id: int,
        status: str,
        progress: float | None = None,
        current_step: str | None = None,
        error_message: str | None = None,
        result_path: str | None = None,
    ) -> Job | None:
        """Update job status and progress."""
        result = await self.session.execute(
            select(Job).where(Job.id == job_id)
        )
        job = result.scalar_one_or_none()
        
        if not job:
            return None
        
        job.status = status
        
        if progress is not None:
            job.progress = progress
        if current_step is not None:
            job.current_step = current_step
        if error_message is not None:
            job.error_message = error_message
        if result_path is not None:
            job.result_path = result_path
        
        # Update timestamps
        now = datetime.now(timezone.utc)
        if status == JobStatus.RUNNING.value and not job.started_at:
            job.started_at = now
        if status in [JobStatus.COMPLETED.value, JobStatus.FAILED.value, JobStatus.CANCELLED.value]:
            job.completed_at = now
        
        await self.session.commit()
        await self.session.refresh(job)
        
        logger.info(f"Job {job_id} status updated: {status}, progress: {progress}")
        return job
    
    async def get_job_with_sample(self, job_id: int) -> tuple[Job, Sample] | None:
        """Get job with its associated sample."""
        result = await self.session.execute(
            select(Job, Sample)
            .join(Sample)
            .where(Job.id == job_id)
        )
        row = result.first()
        return (row[0], row[1]) if row else None
    
    def get_pipeline_config(self, job_type: str) -> dict[str, Any]:
        """Get pipeline configuration for job type."""
        configs = {
            JobType.ALIGNMENT.value: {
                "aligner": "bwa-mem2",
                "threads": settings.pipeline.threads,
                "reference": settings.pipeline.reference_genome,
            },
            JobType.VARIANT_CALLING.value: {
                "caller": "gatk",
                "threads": settings.pipeline.threads,
                "reference": settings.pipeline.reference_genome,
            },
            JobType.ANNOTATION.value: {
                "annotators": ["vep", "clinvar", "gnomad"],
            },
            JobType.QC.value: {
                "tools": ["fastqc", "multiqc"],
            },
            JobType.FULL_PIPELINE.value: {
                "steps": ["qc", "alignment", "variant_calling", "annotation"],
                "threads": settings.pipeline.threads,
                "reference": settings.pipeline.reference_genome,
            },
        }
        return configs.get(job_type, {})
    
    def validate_sample_for_job(self, sample: Sample, job_type: str) -> tuple[bool, str]:
        """Validate that sample has required files for job type."""
        if job_type in [JobType.ALIGNMENT.value, JobType.FULL_PIPELINE.value]:
            if not sample.fastq_r1:
                return False, "Sample missing FASTQ R1 file"
            if not storage_service.file_exists(Path(sample.fastq_r1)):
                return False, "FASTQ R1 file not found on disk"
        
        if job_type == JobType.VARIANT_CALLING.value:
            if not sample.bam_path:
                return False, "Sample missing BAM file (run alignment first)"
            if not storage_service.file_exists(Path(sample.bam_path)):
                return False, "BAM file not found on disk"
        
        if job_type == JobType.ANNOTATION.value:
            if not sample.vcf_path:
                return False, "Sample missing VCF file (run variant calling first)"
            if not storage_service.file_exists(Path(sample.vcf_path)):
                return False, "VCF file not found on disk"
        
        return True, "OK"
