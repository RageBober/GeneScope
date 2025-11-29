"""Project model."""

from datetime import datetime
from typing import TYPE_CHECKING

from sqlalchemy import String, Text, Integer, ForeignKey, DateTime, func
from sqlalchemy.orm import Mapped, mapped_column, relationship

from bioforge.db import Base

if TYPE_CHECKING:
    from .user import User
    from .sample import Sample


class Project(Base):
    """Research project containing samples."""
    
    __tablename__ = "projects"
    
    id: Mapped[int] = mapped_column(primary_key=True)
    name: Mapped[str] = mapped_column(String(255), index=True)
    description: Mapped[str | None] = mapped_column(Text)
    
    # Reference genome
    reference_genome: Mapped[str] = mapped_column(String(50), default="GRCh38")
    
    # Owner
    owner_id: Mapped[int] = mapped_column(ForeignKey("users.id"))
    
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
    owner: Mapped["User"] = relationship(back_populates="projects")
    samples: Mapped[list["Sample"]] = relationship(
        back_populates="project",
        cascade="all, delete-orphan",
    )
    
    def __repr__(self) -> str:
        return f"<Project {self.name}>"
