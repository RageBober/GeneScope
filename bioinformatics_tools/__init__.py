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
"""

__version__ = "0.1.0"
__all__ = [
    "BaseBioTool",
    "ToolConfig",
    "validate_file_path",
    "validate_sequence",
]

from .base_tool import BaseBioTool
from .config import ToolConfig
from .validators import validate_file_path, validate_sequence
