"""
⚡ High-Performance Computing Module for GeneScope

Provides extreme performance optimization for processing hundreds of GB in minutes:

Architecture Layers:
1. **Parallel Processing**: Multiprocessing, Ray distributed computing
2. **GPU Acceleration**: cuDF, Rapids.ai for pandas operations
3. **Streaming**: Chunk-based processing, generator patterns
4. **Distributed**: Dask for big data, Spark-compatible
5. **Caching**: Redis for intermediate results
6. **Optimization**: Numba JIT, vectorization

Performance Targets:
- 100GB FASTQ: <5 minutes (8-core, 64GB RAM)
- 1TB VCF: <15 minutes (16-core, 128GB RAM, GPU)
- 10M variant filtering: <30 seconds (GPU-accelerated)

Security:
- Resource monitoring and limits
- Memory pressure detection
- Automatic cleanup
"""

__version__ = "1.0.0"
__all__ = [
    # Core engines
    "ParallelProcessor",
    "GPUAccelerator",
    "StreamProcessor",
    "DistributedCompute",
    # Utilities
    "PerformanceMonitor",
    "CacheManager",
    "optimize_pipeline",
]

from .parallel import ParallelProcessor
from .gpu import GPUAccelerator
from .streaming import StreamProcessor
from .distributed import DistributedCompute
from .monitoring import PerformanceMonitor
from .cache import CacheManager
from .optimizer import optimize_pipeline
