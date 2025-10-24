"""
🚀 GPU Acceleration Engine

Provides GPU-accelerated data processing using NVIDIA CUDA:
- cuDF for GPU DataFrame operations (100x faster than pandas)
- cuML for GPU machine learning
- RAPIDS ecosystem integration

Performance: 10-100x speedup for large DataFrames (>1M rows)

Requirements:
- NVIDIA GPU with CUDA support
- RAPIDS: pip install cudf-cu11 cuml-cu11
"""

from __future__ import annotations

import logging
from typing import Any, Callable, Optional

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

# Try importing cuDF (GPU-accelerated pandas)
try:
    import cudf
    import cupy as cp
    CUDF_AVAILABLE = True
    logger.info("✅ cuDF available - GPU acceleration enabled")
except ImportError:
    CUDF_AVAILABLE = False
    logger.warning(
        "⚠️  cuDF not available. Install RAPIDS for GPU acceleration:\n"
        "    pip install cudf-cu11 cuml-cu11"
    )

# Try importing cuML (GPU-accelerated scikit-learn)
try:
    import cuml
    CUML_AVAILABLE = True
except ImportError:
    CUML_AVAILABLE = False
    logger.warning("cuML not available for GPU ML")


class GPUAccelerator:
    """
    GPU acceleration for massive dataset processing.

    Automatically falls back to CPU if GPU unavailable.

    Performance Comparison (10M rows):
    - pandas groupby.sum():  ~30 seconds
    - cuDF groupby.sum():    ~0.3 seconds (100x faster)

    Examples:
        >>> accelerator = GPUAccelerator()
        >>> gpu_df = accelerator.to_gpu(df)  # Move to GPU
        >>> result = gpu_df.groupby('category').sum()
        >>> cpu_df = accelerator.to_cpu(result)  # Move back to CPU
    """

    def __init__(self, device_id: int = 0):
        """
        Initialize GPU accelerator.

        Args:
            device_id: GPU device ID (0 for first GPU)
        """
        self.device_id = device_id
        self.gpu_available = CUDF_AVAILABLE

        if self.gpu_available:
            try:
                # Set GPU device
                cp.cuda.Device(device_id).use()
                memory_info = cp.cuda.Device(device_id).mem_info
                total_memory = memory_info[1] / 1e9  # GB
                logger.info(f"GPU {device_id} initialized: {total_memory:.1f} GB memory")
            except Exception as exc:
                logger.error(f"GPU initialization failed: {exc}")
                self.gpu_available = False

    def to_gpu(self, df: pd.DataFrame) -> Any:
        """
        Transfer DataFrame to GPU memory.

        Args:
            df: pandas DataFrame

        Returns:
            cuDF DataFrame (on GPU) or original DataFrame (if no GPU)

        Examples:
            >>> gpu_df = accelerator.to_gpu(df)
        """
        if not self.gpu_available:
            logger.warning("GPU not available, returning CPU DataFrame")
            return df

        try:
            logger.info(f"Transferring DataFrame to GPU ({len(df)} rows, {len(df.columns)} cols)")
            gpu_df = cudf.from_pandas(df)
            logger.info(f"✅ Transfer complete: {self._get_gpu_memory_usage():.2f} GB GPU memory used")
            return gpu_df
        except Exception as exc:
            logger.error(f"GPU transfer failed: {exc}, falling back to CPU")
            return df

    def to_cpu(self, gpu_df: Any) -> pd.DataFrame:
        """
        Transfer DataFrame from GPU to CPU memory.

        Args:
            gpu_df: cuDF DataFrame (on GPU)

        Returns:
            pandas DataFrame

        Examples:
            >>> cpu_df = accelerator.to_cpu(gpu_df)
        """
        if not self.gpu_available or not isinstance(gpu_df, cudf.DataFrame):
            return gpu_df

        try:
            logger.info("Transferring DataFrame from GPU to CPU")
            cpu_df = gpu_df.to_pandas()
            logger.info(f"✅ Transfer complete: {len(cpu_df)} rows")
            return cpu_df
        except Exception as exc:
            logger.error(f"CPU transfer failed: {exc}")
            return gpu_df

    def gpu_groupby_agg(
        self,
        df: pd.DataFrame,
        groupby_cols: list[str],
        agg_dict: dict[str, str | list[str]],
    ) -> pd.DataFrame:
        """
        GPU-accelerated groupby aggregation.

        10-100x faster than pandas for large datasets.

        Args:
            df: Input DataFrame
            groupby_cols: Columns to group by
            agg_dict: Aggregation dictionary {col: ['sum', 'mean', ...]}

        Returns:
            Aggregated DataFrame

        Examples:
            >>> result = accelerator.gpu_groupby_agg(
            ...     df,
            ...     groupby_cols=['category', 'region'],
            ...     agg_dict={'sales': ['sum', 'mean'], 'quantity': 'sum'}
            ... )
        """
        if not self.gpu_available:
            logger.warning("Using CPU groupby (slower)")
            return df.groupby(groupby_cols).agg(agg_dict).reset_index()

        try:
            # Transfer to GPU
            gpu_df = self.to_gpu(df)

            # GPU groupby
            result = gpu_df.groupby(groupby_cols).agg(agg_dict).reset_index()

            # Transfer back to CPU
            return self.to_cpu(result)

        except Exception as exc:
            logger.error(f"GPU groupby failed: {exc}, falling back to CPU")
            return df.groupby(groupby_cols).agg(agg_dict).reset_index()

    def gpu_merge(
        self,
        left: pd.DataFrame,
        right: pd.DataFrame,
        **merge_kwargs: Any,
    ) -> pd.DataFrame:
        """
        GPU-accelerated DataFrame merge/join.

        Much faster for large datasets (>1M rows).

        Args:
            left: Left DataFrame
            right: Right DataFrame
            **merge_kwargs: Arguments for pd.merge()

        Returns:
            Merged DataFrame

        Examples:
            >>> merged = accelerator.gpu_merge(
            ...     df1, df2,
            ...     on='key',
            ...     how='inner'
            ... )
        """
        if not self.gpu_available:
            return pd.merge(left, right, **merge_kwargs)

        try:
            gpu_left = self.to_gpu(left)
            gpu_right = self.to_gpu(right)

            result = gpu_left.merge(gpu_right, **merge_kwargs)

            return self.to_cpu(result)

        except Exception as exc:
            logger.error(f"GPU merge failed: {exc}, falling back to CPU")
            return pd.merge(left, right, **merge_kwargs)

    def gpu_sort_values(
        self,
        df: pd.DataFrame,
        by: str | list[str],
        **sort_kwargs: Any,
    ) -> pd.DataFrame:
        """
        GPU-accelerated sorting.

        Args:
            df: Input DataFrame
            by: Column(s) to sort by
            **sort_kwargs: Arguments for sort_values()

        Returns:
            Sorted DataFrame
        """
        if not self.gpu_available:
            return df.sort_values(by=by, **sort_kwargs)

        try:
            gpu_df = self.to_gpu(df)
            result = gpu_df.sort_values(by=by, **sort_kwargs)
            return self.to_cpu(result)

        except Exception as exc:
            logger.error(f"GPU sort failed: {exc}")
            return df.sort_values(by=by, **sort_kwargs)

    def gpu_apply(
        self,
        df: pd.DataFrame,
        func: Callable,
        **apply_kwargs: Any,
    ) -> pd.DataFrame:
        """
        GPU-accelerated apply (limited support).

        Note: Only vectorized operations benefit from GPU.

        Args:
            df: Input DataFrame
            func: Function to apply
            **apply_kwargs: Arguments for apply()

        Returns:
            DataFrame with function applied
        """
        if not self.gpu_available:
            return df.apply(func, **apply_kwargs)

        try:
            gpu_df = self.to_gpu(df)
            result = gpu_df.apply(func, **apply_kwargs)
            return self.to_cpu(result)

        except Exception as exc:
            logger.error(f"GPU apply failed: {exc}")
            return df.apply(func, **apply_kwargs)

    def _get_gpu_memory_usage(self) -> float:
        """Get current GPU memory usage in GB"""
        if not self.gpu_available:
            return 0.0

        try:
            memory_info = cp.cuda.Device(self.device_id).mem_info
            used = (memory_info[1] - memory_info[0]) / 1e9
            return used
        except Exception:
            return 0.0

    def clear_gpu_memory(self):
        """Clear GPU memory cache"""
        if self.gpu_available:
            try:
                cp.get_default_memory_pool().free_all_blocks()
                logger.info("GPU memory cache cleared")
            except Exception as exc:
                logger.error(f"Failed to clear GPU memory: {exc}")


# ═══════════════════════════════════════════════════════════════
#  GPU-Accelerated Machine Learning
# ═══════════════════════════════════════════════════════════════

class GPUMLAccelerator:
    """
    GPU-accelerated machine learning using cuML.

    Provides drop-in replacement for scikit-learn with 10-50x speedup.

    Examples:
        >>> ml = GPUMLAccelerator()
        >>> model = ml.random_forest_classifier(X_train, y_train)
        >>> predictions = ml.predict(model, X_test)
    """

    def __init__(self):
        """Initialize GPU ML accelerator"""
        self.gpu_available = CUML_AVAILABLE

        if not self.gpu_available:
            logger.warning("cuML not available, using CPU scikit-learn")

    def random_forest_classifier(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        **rf_kwargs: Any,
    ) -> Any:
        """
        GPU-accelerated Random Forest classifier.

        Args:
            X: Feature matrix
            y: Target labels
            **rf_kwargs: RandomForestClassifier parameters

        Returns:
            Trained model
        """
        if self.gpu_available:
            from cuml.ensemble import RandomForestClassifier as cuRF
            model = cuRF(**rf_kwargs)
        else:
            from sklearn.ensemble import RandomForestClassifier
            model = RandomForestClassifier(**rf_kwargs)

        model.fit(X, y)
        return model

    def predict(self, model: Any, X: pd.DataFrame) -> np.ndarray:
        """
        Predict using trained model.

        Args:
            model: Trained model
            X: Feature matrix

        Returns:
            Predictions
        """
        return model.predict(X)


# ═══════════════════════════════════════════════════════════════
#  Convenience Functions
# ═══════════════════════════════════════════════════════════════

def gpu_groupby(
    df: pd.DataFrame,
    groupby_cols: list[str],
    agg_dict: dict[str, Any],
) -> pd.DataFrame:
    """
    Quick GPU groupby (convenience function).

    Examples:
        >>> result = gpu_groupby(
        ...     df,
        ...     ['category'],
        ...     {'sales': 'sum'}
        ... )
    """
    accelerator = GPUAccelerator()
    return accelerator.gpu_groupby_agg(df, groupby_cols, agg_dict)


def gpu_merge(
    left: pd.DataFrame,
    right: pd.DataFrame,
    **merge_kwargs: Any,
) -> pd.DataFrame:
    """Quick GPU merge (convenience function)"""
    accelerator = GPUAccelerator()
    return accelerator.gpu_merge(left, right, **merge_kwargs)
