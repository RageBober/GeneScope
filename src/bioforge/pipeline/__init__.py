"""BioForge Pipeline Module.

Genomic analysis pipeline components:
- Alignment (BWA, Minimap2, STAR)
- Variant Calling (GATK, DeepVariant)
- Annotation (VEP, ClinVar, gnomAD)
- Quality Control (FastQC, MultiQC)
"""

from .alignment import AlignmentEngine, AlignmentStats
from .variant_calling import VariantCaller, VariantStats, JointGenotyping
from .annotation import VariantAnnotator, AnnotatedVariant
from .qc import QualityController, QCMetrics
from .validators import (
    FileValidator,
    ParameterValidator,
    OutputValidator,
    PipelineValidator,
)

__all__ = [
    # Alignment
    "AlignmentEngine",
    "AlignmentStats",
    # Variant Calling
    "VariantCaller",
    "VariantStats",
    "JointGenotyping",
    # Annotation
    "VariantAnnotator",
    "AnnotatedVariant",
    # QC
    "QualityController",
    "QCMetrics",
    # Validators
    "FileValidator",
    "ParameterValidator",
    "OutputValidator",
    "PipelineValidator",
]
