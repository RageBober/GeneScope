"""Parallel processing module for genomic data analysis.

This module provides distributed and parallel processing capabilities
for large-scale genomic data using Dask and other parallel frameworks.
"""

from .chunk_managers import BaseChunkManager
from .chunk_managers import CSVChunkManager
from .chunk_managers import VCFChunkManager
from .chunk_managers import get_chunk_manager
from .dask_processor import DaskGenomicProcessor
from .performance_monitor import PerformanceMonitor
from .performance_monitor import TaskMetrics

__all__ = [
    "BaseChunkManager",
    "CSVChunkManager",
    "DaskGenomicProcessor",
    "PerformanceMonitor",
    "TaskMetrics",
    "VCFChunkManager",
    "get_chunk_manager",
]
