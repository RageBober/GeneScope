"""Sample model."""

from datetime import datetime
from enum import Enum
from typing import TYPE_CHECKING

from sqlalchemy import String, Text, Integer, ForeignKey, DateTime, JSON, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from bioforge.db import Base

if TYPE_CHECKING:
    from .project import Project
    from .job import Job
    from .variant import Variant


class SampleStatus(str, Enum):
    """Sample processing status."""
    PENDING = "pending"
    UPLOADING = "uploading"
    UPLOADED = "uploaded"
    PROCESSING = "processing"
    COMPLETED = "completed"
    FAILED = "failed"


class Sample(Base):
    """Biological sample with sequencing data."""
    
    __tablename__ = "samples"
    
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(255), index=True)
    description: Mapped[str | None] = mapped_column(Text)
    
    # Sample metadata
    sample_type: Mapped[str | None] = mapped_column(String(50))  # WGS, WES, RNA-seq
    organism: Mapped[str] = mapped_column(String(100), default="Homo sapiens")
    tissue: Mapped[str | None] = mapped_column(String(100))
    
    # File paths (relative to storage)
    fastq_r1: Mapped[str | None] = mapped_column(String(500))
    fastq_r2: Mapped[str | None] = mapped_column(String(500))
    bam_path: Mapped[str | None] = mapped_column(String(500))
    vcf_path: Mapped[str | None] = mapped_column(String(500))
    
    # Status
    status: Mapped[str] = mapped_column(String(50), default=SampleStatus.PENDING.value)
    
    # QC metrics (JSON)
    qc_metrics: Mapped[dict | None] = mapped_column(JSON)
    
    # Project
    project_id: Mapped[int] = mapped_column(ForeignKey("projects.id"))
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        server_default=func.now(),
        onupdate=func.now(),
    )
    
    # Relationships
    project: Mapped["Project"] = relationship(back_populates="samples")
    jobs: Mapped[list["Job"]] = relationship(
        back_populates="sample",
        cascade="all, delete-orphan",
    )
    variants: Mapped[list["Variant"]] = relationship(
        back_populates="sample",
        cascade="all, delete-orphan",
    )
    
    def __repr__(self) -> str:
        return f"<Sample {self.name}>"
