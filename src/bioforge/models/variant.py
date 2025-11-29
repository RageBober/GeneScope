"""Variant model."""

from datetime import datetime, timezone
from typing import TYPE_CHECKING

from sqlalchemy import DateTime, Float, ForeignKey, Index, Integer, String, Text, JSON
from sqlalchemy.orm import Mapped, mapped_column, relationship

from bioforge.db import Base

if TYPE_CHECKING:
    from bioforge.models import Sample


class Variant(Base):
    """
    Genetic variant detected in a sample.
    
    Stores SNPs, indels, and structural variants with annotations.
    """
    
    __tablename__ = "variants"
    
    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    
    # Location
    chromosome: Mapped[str] = mapped_column(String(10), nullable=False)
    position: Mapped[int] = mapped_column(Integer, nullable=False)
    ref: Mapped[str] = mapped_column(String(500), nullable=False)
    alt: Mapped[str] = mapped_column(String(500), nullable=False)
    
    # Variant info
    variant_type: Mapped[str] = mapped_column(String(20), default="SNV")  # SNV, INDEL, SV
    quality: Mapped[float | None] = mapped_column(Float, nullable=True)
    filter_status: Mapped[str] = mapped_column(String(50), default="PASS")
    
    # Genotype
    genotype: Mapped[str | None] = mapped_column(String(10), nullable=True)  # 0/1, 1/1, etc.
    depth: Mapped[int | None] = mapped_column(Integer, nullable=True)
    allele_depth: Mapped[str | None] = mapped_column(String(50), nullable=True)  # e.g., "30,15"
    
    # Gene info
    gene_symbol: Mapped[str | None] = mapped_column(String(50), nullable=True, index=True)
    gene_id: Mapped[str | None] = mapped_column(String(50), nullable=True)
    transcript_id: Mapped[str | None] = mapped_column(String(50), nullable=True)
    
    # Effect
    consequence: Mapped[str | None] = mapped_column(String(100), nullable=True)
    impact: Mapped[str | None] = mapped_column(String(20), nullable=True)  # HIGH, MODERATE, LOW, MODIFIER
    hgvs_c: Mapped[str | None] = mapped_column(String(200), nullable=True)  # c.1234A>G
    hgvs_p: Mapped[str | None] = mapped_column(String(200), nullable=True)  # p.Arg123Gly
    
    # Population frequencies
    gnomad_af: Mapped[float | None] = mapped_column(Float, nullable=True)
    gnomad_af_popmax: Mapped[float | None] = mapped_column(Float, nullable=True)
    
    # Clinical
    clinvar_id: Mapped[str | None] = mapped_column(String(20), nullable=True)
    clinvar_significance: Mapped[str | None] = mapped_column(String(100), nullable=True)
    clinvar_conditions: Mapped[str | None] = mapped_column(Text, nullable=True)
    
    # Pathogenicity scores
    cadd_score: Mapped[float | None] = mapped_column(Float, nullable=True)
    revel_score: Mapped[float | None] = mapped_column(Float, nullable=True)
    sift_pred: Mapped[str | None] = mapped_column(String(20), nullable=True)
    polyphen_pred: Mapped[str | None] = mapped_column(String(20), nullable=True)
    
    # Additional annotations (flexible JSON)
    annotations: Mapped[dict | None] = mapped_column(JSON, nullable=True)
    
    # Relationships
    sample_id: Mapped[int] = mapped_column(ForeignKey("samples.id", ondelete="CASCADE"))
    sample: Mapped["Sample"] = relationship(back_populates="variants")
    
    # Timestamps
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True),
        default=lambda: datetime.now(timezone.utc),
    )
    
    # Indexes for fast queries
    __table_args__ = (
        Index("ix_variants_location", "chromosome", "position"),
        Index("ix_variants_clinvar", "clinvar_significance"),
        Index("ix_variants_impact", "impact"),
        Index("ix_variants_sample_gene", "sample_id", "gene_symbol"),
    )
    
    def __repr__(self) -> str:
        return f"<Variant {self.chromosome}:{self.position} {self.ref}>{self.alt}>"
    
    @property
    def variant_id(self) -> str:
        """Standard variant ID format."""
        return f"{self.chromosome}-{self.position}-{self.ref}-{self.alt}"
    
    @property
    def is_pathogenic(self) -> bool:
        """Check if variant is pathogenic based on ClinVar."""
        if not self.clinvar_significance:
            return False
        return "pathogenic" in self.clinvar_significance.lower()
    
    @property
    def is_rare(self, threshold: float = 0.01) -> bool:
        """Check if variant is rare (AF < 1%)."""
        if self.gnomad_af is None:
            return True  # Unknown = potentially rare
        return self.gnomad_af < threshold
