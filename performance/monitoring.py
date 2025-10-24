"""
📊 Performance Monitoring

Track and optimize pipeline performance:
- Execution time tracking
- Memory usage monitoring
- Throughput calculation
- Performance profiling
"""

from __future__ import annotations

import logging
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Optional

import psutil

logger = logging.getLogger(__name__)


@dataclass
class PerformanceMetrics:
    """Performance metrics for an operation"""
    operation_name: str
    duration_seconds: float
    memory_used_mb: float
    throughput_mb_per_sec: float = 0.0
    rows_processed: int = 0
    rows_per_second: float = 0.0
    metadata: dict[str, Any] = field(default_factory=dict)


class PerformanceMonitor:
    """
    Monitor and track performance metrics.

    Examples:
        >>> monitor = PerformanceMonitor()
        >>> with monitor.track("processing"):
        ...     # Do expensive operation
        ...     process_data()
        >>> metrics = monitor.get_metrics("processing")
        >>> print(f"Duration: {metrics.duration_seconds:.2f}s")
    """

    def __init__(self):
        """Initialize performance monitor"""
        self.metrics: dict[str, PerformanceMetrics] = {}
        self.process = psutil.Process()

    @contextmanager
    def track(
        self,
        operation_name: str,
        rows: int = 0,
        data_size_mb: float = 0.0,
    ):
        """
        Context manager to track operation performance.

        Args:
            operation_name: Name of operation
            rows: Number of rows processed
            data_size_mb: Size of data in MB

        Examples:
            >>> with monitor.track("filtering", rows=1000000):
            ...     filtered = filter_data(df)
        """
        logger.info(f"🚀 Starting: {operation_name}")

        # Record start state
        start_time = time.time()
        start_memory = self.process.memory_info().rss / 1024 / 1024  # MB

        try:
            yield
        finally:
            # Record end state
            end_time = time.time()
            end_memory = self.process.memory_info().rss / 1024 / 1024  # MB

            duration = end_time - start_time
            memory_used = end_memory - start_memory

            # Calculate throughput
            throughput = data_size_mb / duration if duration > 0 and data_size_mb > 0 else 0.0
            rows_per_sec = rows / duration if duration > 0 and rows > 0 else 0.0

            # Store metrics
            self.metrics[operation_name] = PerformanceMetrics(
                operation_name=operation_name,
                duration_seconds=duration,
                memory_used_mb=memory_used,
                throughput_mb_per_sec=throughput,
                rows_processed=rows,
                rows_per_second=rows_per_sec,
            )

            logger.info(
                f"✅ Completed: {operation_name} "
                f"({duration:.2f}s, {memory_used:.1f}MB, "
                f"{rows_per_sec:.0f} rows/s)"
            )

    def get_metrics(self, operation_name: str) -> Optional[PerformanceMetrics]:
        """Get metrics for an operation"""
        return self.metrics.get(operation_name)

    def get_all_metrics(self) -> dict[str, PerformanceMetrics]:
        """Get all metrics"""
        return self.metrics

    def print_summary(self):
        """Print performance summary"""
        print("\n" + "=" * 60)
        print("PERFORMANCE SUMMARY")
        print("=" * 60)

        for name, metrics in self.metrics.items():
            print(f"\n{name}:")
            print(f"  Duration:    {metrics.duration_seconds:.2f}s")
            print(f"  Memory:      {metrics.memory_used_mb:.1f} MB")
            if metrics.rows_processed > 0:
                print(f"  Throughput:  {metrics.rows_per_second:.0f} rows/s")
            if metrics.throughput_mb_per_sec > 0:
                print(f"  Data rate:   {metrics.throughput_mb_per_sec:.1f} MB/s")

        print("=" * 60 + "\n")
