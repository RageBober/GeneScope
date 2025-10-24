"""
⚡ Pipeline Optimizer

Automatic performance optimization:
- Data size detection
- Strategy selection (CPU/GPU/Distributed)
- Memory optimization
- Numba JIT compilation
"""

from __future__ import annotations

import logging
from typing import Any, Callable

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

# Try importing Numba
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(*args, **kwargs):
        """Dummy decorator if Numba not available"""
        def decorator(func):
            return func
        return decorator if not args else decorator(args[0])


def optimize_pipeline(
    data_size_gb: float,
    n_workers: int = 8,
    gpu_available: bool = False,
) -> dict[str, Any]:
    """
    Recommend optimal processing strategy based on data size.

    Args:
        data_size_gb: Dataset size in GB
        n_workers: Available CPU cores
        gpu_available: Whether GPU is available

    Returns:
        Dictionary with recommendations

    Examples:
        >>> recommendations = optimize_pipeline(
        ...     data_size_gb=50,
        ...     n_workers=16,
        ...     gpu_available=True
        ... )
        >>> print(recommendations['strategy'])  # 'gpu' or 'distributed'
    """
    recommendations = {
        'data_size_gb': data_size_gb,
        'strategy': None,
        'chunk_size': None,
        'use_parallel': False,
        'use_gpu': False,
        'use_distributed': False,
        'estimated_time_minutes': None,
    }

    # Small data (<1GB): Single process
    if data_size_gb < 1:
        recommendations['strategy'] = 'single_process'
        recommendations['chunk_size'] = 100_000
        recommendations['estimated_time_minutes'] = 1

    # Medium data (1-10GB): Multiprocessing
    elif data_size_gb < 10:
        recommendations['strategy'] = 'multiprocessing'
        recommendations['use_parallel'] = True
        recommendations['chunk_size'] = 50_000
        recommendations['estimated_time_minutes'] = data_size_gb * 2

    # Large data with GPU (10-100GB): GPU acceleration
    elif data_size_gb < 100 and gpu_available:
        recommendations['strategy'] = 'gpu'
        recommendations['use_gpu'] = True
        recommendations['chunk_size'] = 1_000_000
        recommendations['estimated_time_minutes'] = data_size_gb * 0.5

    # Large data without GPU (10-100GB): Distributed
    elif data_size_gb < 100:
        recommendations['strategy'] = 'distributed'
        recommendations['use_distributed'] = True
        recommendations['chunk_size'] = 100_000
        recommendations['estimated_time_minutes'] = data_size_gb * 3

    # Huge data (>100GB): Distributed + Streaming
    else:
        recommendations['strategy'] = 'distributed_streaming'
        recommendations['use_distributed'] = True
        recommendations['chunk_size'] = 50_000
        recommendations['estimated_time_minutes'] = data_size_gb * 2

    logger.info(
        f"📊 Optimization: {data_size_gb:.1f}GB → {recommendations['strategy']} "
        f"(~{recommendations['estimated_time_minutes']:.0f} min)"
    )

    return recommendations


def optimize_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Optimize DataFrame memory usage.

    Converts dtypes to more efficient representations.

    Args:
        df: Input DataFrame

    Returns:
        Optimized DataFrame

    Examples:
        >>> df_optimized = optimize_dataframe(df)
        >>> print(f"Memory saved: {df.memory_usage().sum() - df_optimized.memory_usage().sum()}")
    """
    initial_memory = df.memory_usage(deep=True).sum() / 1024 / 1024  # MB

    # Optimize numeric columns
    for col in df.select_dtypes(include=['int']).columns:
        col_min = df[col].min()
        col_max = df[col].max()

        if col_min >= 0:
            if col_max < 256:
                df[col] = df[col].astype(np.uint8)
            elif col_max < 65536:
                df[col] = df[col].astype(np.uint16)
            elif col_max < 4294967296:
                df[col] = df[col].astype(np.uint32)
        else:
            if col_min > -128 and col_max < 128:
                df[col] = df[col].astype(np.int8)
            elif col_min > -32768 and col_max < 32768:
                df[col] = df[col].astype(np.int16)
            elif col_min > -2147483648 and col_max < 2147483648:
                df[col] = df[col].astype(np.int32)

    # Optimize float columns
    for col in df.select_dtypes(include=['float']).columns:
        df[col] = df[col].astype(np.float32)

    # Convert object to category if beneficial
    for col in df.select_dtypes(include=['object']).columns:
        num_unique = df[col].nunique()
        num_total = len(df[col])

        if num_unique / num_total < 0.5:  # Less than 50% unique
            df[col] = df[col].astype('category')

    final_memory = df.memory_usage(deep=True).sum() / 1024 / 1024  # MB
    reduction = (1 - final_memory / initial_memory) * 100

    logger.info(
        f"💾 Memory optimized: {initial_memory:.1f}MB → {final_memory:.1f}MB "
        f"({reduction:.1f}% reduction)"
    )

    return df


@jit(nopython=True)
def fast_mean(arr: np.ndarray) -> float:
    """Numba-accelerated mean calculation"""
    return np.mean(arr)


@jit(nopython=True)
def fast_sum(arr: np.ndarray) -> float:
    """Numba-accelerated sum calculation"""
    return np.sum(arr)
