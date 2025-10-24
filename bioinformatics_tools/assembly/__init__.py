"""
🧩 Genome Assembly Tools Module

Provides wrappers for:
- MEGAHIT: Ultra-fast and memory-efficient NGS assembler
- SPAdes: Versatile genome assembler for single-cell and standard assemblies
"""

__all__ = ["MEGAHITTool", "SPAdesTool"]

from .megahit import MEGAHITTool
from .spades import SPAdesTool
