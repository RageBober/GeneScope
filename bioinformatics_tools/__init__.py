"""
🧬 Bioinformatics Tools Integration Module

Provides secure, containerized wrappers for external bioinformatics tools:
- MetaGraph: DNA/RNA sequence indexing and search
- Kraken2/Centrifuge: Taxonomic classification
- MEGAHIT/SPAdes: Genome assembly
- GTDB-Tk: Bacterial classification
- XtriMoPGLM: Protein property prediction
- DeepSEA/Enformer: Mutation effect prediction

Security features:
- Input validation and sanitization
- Docker containerization
- Resource limiting
- Path traversal prevention

BioForge Integration:
- MetagenomicsPipeline: Kraken2 → Assembly → GTDB-Tk
- ProteinAnalysisPipeline: XtriMoPGLM protein prediction
- MutationEffectPipeline: Enformer mutation effects
- MetaGraphDataSource: Global sequence search
"""

__version__ = "0.2.0"
__all__ = [
    # Core components
    "BaseBioTool",
    "ToolConfig",
    "validate_file_path",
    "validate_sequence",
    # BioForge pipelines
    "MetagenomicsPipeline",
    "ProteinAnalysisPipeline",
    "MutationEffectPipeline",
    "MetaGraphDataSource",
    "PipelineResult",
]

from .base_tool import BaseBioTool
from .config import ToolConfig
from .validators import validate_file_path, validate_sequence
from .bioforge_integration import (
    MetagenomicsPipeline,
    ProteinAnalysisPipeline,
    MutationEffectPipeline,
    MetaGraphDataSource,
    PipelineResult,
)
