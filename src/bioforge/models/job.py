"""Job model for pipeline execution tracking."""

from datetime import datetime
from enum import Enum
from typing import TYPE_CHECKING

from sqlalchemy import String, Text, Integer, Float, ForeignKey, DateTime, JSON, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from bioforge.db import Base

if TYPE_CHECKING:
    from .sample import Sample


class JobStatus(str, Enum):
    """Job execution status."""
    PENDING = "pending"
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobType(str, Enum):
    """Type of pipeline job."""
    ALIGNMENT = "alignment"
    VARIANT_CALLING = "variant_calling"
    ANNOTATION = "annotation"
    QC = "qc"
    FULL_PIPELINE = "full_pipeline"


class Job(Base):
    """Pipeline job execution."""
    
    __tablename__ = "jobs"
    
    id: Mapped[int] = mapped_column(primary_key=True)
    
    # Job type and status
    job_type: Mapped[str] = mapped_column(String(50))
    status: Mapped[str] = mapped_column(String(50), default=JobStatus.PENDING.value, index=True)
    
    # Progress
    progress: Mapped[float] = mapped_column(Float, default=0.0)  # 0-100
    current_step: Mapped[str | None] = mapped_column(String(100))
    
    # Results
    result_path: Mapped[str | None] = mapped_column(String(500))
    result_data: Mapped[dict | None] = mapped_column(JSON)  # Structured results (stats, paths, etc.)
    error_message: Mapped[str | None] = mapped_column(Text)
    
    # Parameters and logs (JSON)
    parameters: Mapped[dict | None] = mapped_column(JSON)
    logs: Mapped[list | None] = mapped_column(JSON)
    
    # Task queue ID (for ARQ)
    task_id: Mapped[str | None] = mapped_column(String(100), index=True)
    
    # Sample
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id"))
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
    )
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    completed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    
    # Relationships
    sample: Mapped["Sample"] = relationship(back_populates="jobs")
    
    @property
    def duration_seconds(self) -> float | None:
        """Calculate job duration in seconds."""
        if self.started_at and self.completed_at:
            return (self.completed_at - self.started_at).total_seconds()
        return None
    
    def __repr__(self) -> str:
        return f"<Job {self.id} {self.job_type} {self.status}>"
