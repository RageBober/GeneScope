"""Dask-based parallel processor for genomic data analysis.

This module provides distributed processing capabilities for large genomic datasets
using Dask for parallel computation and memory-efficient processing.
"""

import logging
import time
from pathlib import Path
from typing import Any

import pandas as pd
from dask import delayed
from dask.distributed import Client
from dask.distributed import as_completed
from dask.distributed import get_client

from ..core.logging_config import get_logger
from .chunk_managers import get_chunk_manager
from .performance_monitor import PerformanceMonitor

logger = get_logger(__name__)


class DaskGenomicProcessor:
    """Parallel genomic data processor using Dask for distributed computing.
    
    This class provides high-performance processing of large genomic datasets
    through intelligent chunking and distributed computation.
    """

    def __init__(self,
                 n_workers: int = 4,
                 memory_limit: str = "2GB",
                 threads_per_worker: int = 2,
                 dashboard_address: str | None = ":8787"):
        """Initialize Dask genomic processor.
        
        Args:
            n_workers: Number of worker processes
            memory_limit: Memory limit per worker
            threads_per_worker: Number of threads per worker
            dashboard_address: Dask dashboard address (None to disable)
        """
        self.n_workers = n_workers
        self.memory_limit = memory_limit
        self.threads_per_worker = threads_per_worker
        self.dashboard_address = dashboard_address

        self.client: Client | None = None
        self.performance_monitor = PerformanceMonitor()

        # Chunk managers will be created on demand
        self._chunk_managers = {}

        self._initialize_cluster()

    def _initialize_cluster(self) -> None:
        """Initialize Dask distributed cluster."""
        try:
            # Try to connect to existing cluster first
            self.client = get_client()
            logger.info("Connected to existing Dask cluster")
        except (ImportError, ValueError):
            # Create new local cluster
            try:
                self.client = Client(
                    n_workers=self.n_workers,
                    threads_per_worker=self.threads_per_worker,
                    memory_limit=self.memory_limit,
                    dashboard_address=self.dashboard_address,
                    silence_logs=logging.WARNING
                )
                logger.info(f"Created Dask cluster with {self.n_workers} workers")
                logger.info(f"Dashboard available at: {self.client.dashboard_link}")
            except Exception as e:
                logger.error(f"Failed to initialize Dask cluster: {e}")
                self.client = None

    @property
    def is_distributed(self) -> bool:
        """Check if distributed processing is available."""
        return self.client is not None

    @delayed
    def process_genomic_chunk(self,
                            chunk_data: pd.DataFrame,
                            analysis_type: str = "basic",
                            **kwargs) -> dict[str, Any]:
        """Process a single genomic data chunk.
        
        Args:
            chunk_data: Chunk of genomic data
            analysis_type: Type of analysis to perform
            **kwargs: Additional analysis parameters
            
        Returns:
            Dictionary with analysis results
        """
        start_time = time.time()

        try:
            if analysis_type == "variant_stats":
                results = self._compute_variant_statistics(chunk_data, **kwargs)
            elif analysis_type == "quality_metrics":
                results = self._compute_quality_metrics(chunk_data, **kwargs)
            elif analysis_type == "annotation":
                results = self._annotate_variants(chunk_data, **kwargs)
            elif analysis_type == "filtering":
                results = self._filter_variants(chunk_data, **kwargs)
            else:
                results = self._basic_analysis(chunk_data, **kwargs)

            processing_time = time.time() - start_time
            results["processing_time"] = processing_time
            results["chunk_size"] = len(chunk_data)

            return results

        except Exception as e:
            logger.error(f"Error processing chunk: {e}")
            return {
                "error": str(e),
                "chunk_size": len(chunk_data) if chunk_data is not None else 0,
                "processing_time": time.time() - start_time
            }

    def process_large_file_parallel(self,
                                  file_path: str | Path,
                                  file_type: str = "csv",
                                  analysis_type: str = "basic",
                                  chunk_size_mb: int = 100,
                                  **analysis_kwargs) -> dict[str, Any]:
        """Process large genomic file in parallel chunks.
        
        Args:
            file_path: Path to input file
            file_type: Type of file (csv, vcf, etc.)
            analysis_type: Type of analysis to perform
            chunk_size_mb: Size of chunks in MB
            **analysis_kwargs: Analysis parameters
            
        Returns:
            Aggregated results from all chunks
        """
        file_path = Path(file_path)

        if not file_path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")

        file_size_mb = file_path.stat().st_size / (1024 * 1024)
        logger.info(f"Processing {file_size_mb:.2f}MB file with {self.n_workers} workers")

        # Start performance monitoring
        monitor_task = self.performance_monitor.start_monitoring(
            task_name=f"process_{file_path.name}",
            expected_duration=file_size_mb * 2  # rough estimate: 2 sec per MB
        )

        try:
            # Create chunks using appropriate chunk manager
            chunk_manager = get_chunk_manager(file_type)
            chunks = chunk_manager.create_chunks(file_path, chunk_size_mb)

            logger.info(f"Created {len(chunks)} chunks for processing")

            if not self.is_distributed:
                logger.warning("Dask cluster not available, falling back to sequential processing")
                return self._process_sequential(chunks, analysis_type, **analysis_kwargs)

            # Process chunks in parallel
            futures = []
            for i, chunk in enumerate(chunks):
                future = self.process_genomic_chunk(
                    chunk, analysis_type, chunk_id=i, **analysis_kwargs
                )
                futures.append(future)

            # Collect results
            results = []
            completed_futures = as_completed(futures)

            for future in completed_futures:
                try:
                    result = future.result()
                    results.append(result)
                    logger.debug(f"Completed chunk processing: {len(results)}/{len(chunks)}")
                except Exception as e:
                    logger.error(f"Failed to process chunk: {e}")
                    results.append({"error": str(e)})

            # Aggregate results
            aggregated_results = self._aggregate_results(results, analysis_type)

            # Stop monitoring
            self.performance_monitor.stop_monitoring(monitor_task)

            # Add performance metrics
            performance_stats = self.performance_monitor.get_task_stats(monitor_task)
            aggregated_results["performance"] = performance_stats

            return aggregated_results

        except Exception as e:
            self.performance_monitor.stop_monitoring(monitor_task)
            logger.error(f"Error in parallel processing: {e}")
            raise

    def _basic_analysis(self, data: pd.DataFrame, **kwargs) -> dict[str, Any]:
        """Basic analysis for any genomic data chunk."""
        try:
            analysis = {
                "record_count": len(data),
                "columns": list(data.columns),
                "data_types": data.dtypes.to_dict(),
                "memory_usage_mb": data.memory_usage(deep=True).sum() / (1024 * 1024)
            }

            # Basic statistics for numeric columns
            numeric_cols = data.select_dtypes(include=["number"]).columns
            if len(numeric_cols) > 0:
                analysis["numeric_summary"] = data[numeric_cols].describe().to_dict()

            return analysis

        except Exception as e:
            logger.error(f"Error in basic analysis: {e}")
            return {"error": str(e)}

    def _compute_variant_statistics(self, data: pd.DataFrame, **kwargs) -> dict[str, Any]:
        """Compute variant statistics for a chunk."""
        try:
            stats = {}

            if "CHROM" in data.columns:
                stats["chromosomes"] = data["CHROM"].value_counts().to_dict()

            if "QUAL" in data.columns:
                stats["quality_stats"] = {
                    "mean": float(data["QUAL"].mean()),
                    "median": float(data["QUAL"].median()),
                    "std": float(data["QUAL"].std())
                }

            if "REF" in data.columns and "ALT" in data.columns:
                # Variant type analysis
                data["variant_type"] = data.apply(
                    lambda row: self._classify_variant_type(row["REF"], str(row["ALT"])),
                    axis=1
                )
                stats["variant_types"] = data["variant_type"].value_counts().to_dict()

            stats["total_variants"] = len(data)

            return stats

        except Exception as e:
            logger.error(f"Error computing variant statistics: {e}")
            return {"error": str(e)}

    def _compute_quality_metrics(self, data: pd.DataFrame, **kwargs) -> dict[str, Any]:
        """Compute quality metrics for genomic data chunk."""
        try:
            metrics = {}

            # Basic data quality
            metrics["total_records"] = len(data)
            metrics["missing_values"] = data.isnull().sum().to_dict()
            metrics["duplicate_records"] = int(data.duplicated().sum())

            # Genomic-specific quality metrics
            if "QUAL" in data.columns:
                qual_data = pd.to_numeric(data["QUAL"], errors="coerce")
                metrics["quality_distribution"] = {  # type: ignore[assignment]
                    "high_quality": int((qual_data >= 30).sum()),
                    "medium_quality": int(((qual_data >= 10) & (qual_data < 30)).sum()),
                    "low_quality": int((qual_data < 10).sum())
                }

            if "DP" in data.columns:  # Depth of coverage
                dp_data = pd.to_numeric(data["DP"], errors="coerce")
                metrics["coverage_stats"] = {  # type: ignore[assignment]
                    "mean_depth": float(dp_data.mean()) if not dp_data.empty else 0,
                    "low_coverage": int((dp_data < 10).sum()),
                    "high_coverage": int((dp_data >= 30).sum())
                }

            return metrics

        except Exception as e:
            logger.error(f"Error computing quality metrics: {e}")
            return {"error": str(e)}

    def _annotate_variants(self, data: pd.DataFrame, **kwargs) -> dict[str, Any]:
        """Annotate variants in a chunk (placeholder for future enhancement)."""
        try:
            # Basic annotation - in the future integrate with external databases
            annotations = {}

            if "CHROM" in data.columns and "POS" in data.columns:
                # Create basic genomic coordinates annotation
                data["genomic_coordinate"] = data.apply(
                    lambda row: f"{row['CHROM']}:{row['POS']}", axis=1
                )
                annotations["annotated_count"] = len(data)

            # Placeholder for more sophisticated annotation
            annotations["annotation_source"] = "basic"  # type: ignore[assignment]

            return annotations

        except Exception as e:
            logger.error(f"Error annotating variants: {e}")
            return {"error": str(e)}

    def _filter_variants(self, data: pd.DataFrame, quality_threshold: float = 20.0, **kwargs) -> dict[str, Any]:
        """Filter variants based on quality criteria."""
        try:
            original_count = len(data)

            # Apply quality filter
            if "QUAL" in data.columns:
                qual_data = pd.to_numeric(data["QUAL"], errors="coerce")
                filtered_data = data[qual_data >= quality_threshold]
            else:
                filtered_data = data

            filtered_count = len(filtered_data)

            return {
                "original_count": original_count,
                "filtered_count": filtered_count,
                "removed_count": original_count - filtered_count,
                "filter_criteria": f"QUAL >= {quality_threshold}"
            }

        except Exception as e:
            logger.error(f"Error filtering variants: {e}")
            return {"error": str(e)}

    def _classify_variant_type(self, ref: str, alt: str) -> str:
        """Classify variant type based on REF and ALT alleles."""
        try:
            if len(ref) == len(alt) == 1:
                return "SNV"  # Single Nucleotide Variant
            if len(ref) > len(alt):
                return "DEL"  # Deletion
            if len(ref) < len(alt):
                return "INS"  # Insertion
            return "COMPLEX"  # Complex variant
        except Exception:
            return "UNKNOWN"

    def _aggregate_results(self, results: list[dict[str, Any]], analysis_type: str) -> dict[str, Any]:
        """Aggregate results from multiple chunks."""
        try:
            aggregated = {
                "total_chunks": len(results),
                "successful_chunks": len([r for r in results if "error" not in r]),
                "failed_chunks": len([r for r in results if "error" in r]),
                "analysis_type": analysis_type
            }

            # Filter out failed chunks for aggregation
            successful_results = [r for r in results if "error" not in r]

            if not successful_results:
                aggregated["error"] = "No chunks processed successfully"
                return aggregated

            # Aggregate based on analysis type
            if analysis_type == "variant_stats":
                aggregated.update(self._aggregate_variant_stats(successful_results))
            elif analysis_type == "quality_metrics":
                aggregated.update(self._aggregate_quality_metrics(successful_results))
            elif analysis_type == "filtering":
                aggregated.update(self._aggregate_filtering_results(successful_results))
            else:
                aggregated.update(self._aggregate_basic_results(successful_results))

            # Performance aggregation
            total_processing_time = sum(r.get("processing_time", 0) for r in results)
            total_records = sum(r.get("chunk_size", 0) for r in results)

            aggregated["performance_summary"] = {
                "total_processing_time_seconds": total_processing_time,
                "total_records_processed": total_records,
                "average_processing_speed": total_records / total_processing_time if total_processing_time > 0 else 0,
                "parallel_efficiency": (total_processing_time / len(results)) if len(results) > 0 else 0
            }

            return aggregated

        except Exception as e:
            logger.error(f"Error aggregating results: {e}")
            return {
                "error": f"Aggregation failed: {e!s}",
                "raw_results": results
            }

    def _aggregate_variant_stats(self, results: list[dict[str, Any]]) -> dict[str, Any]:
        """Aggregate variant statistics from multiple chunks."""
        aggregated_stats = {}

        # Aggregate chromosome counts
        chromosome_counts = {}
        variant_type_counts = {}
        total_variants = 0

        quality_stats_list = []

        for result in results:
            total_variants += result.get("total_variants", 0)

            # Chromosome counts
            if "chromosomes" in result:
                for chrom, count in result["chromosomes"].items():
                    chromosome_counts[chrom] = chromosome_counts.get(chrom, 0) + count

            # Variant types
            if "variant_types" in result:
                for vtype, count in result["variant_types"].items():
                    variant_type_counts[vtype] = variant_type_counts.get(vtype, 0) + count

            # Quality stats
            if "quality_stats" in result:
                quality_stats_list.append(result["quality_stats"])

        aggregated_stats["total_variants"] = total_variants  # type: ignore[assignment]
        aggregated_stats["chromosomes"] = chromosome_counts  # type: ignore[assignment]
        aggregated_stats["variant_types"] = variant_type_counts  # type: ignore[assignment]

        # Aggregate quality statistics
        if quality_stats_list:
            aggregated_stats["quality_stats"] = {  # type: ignore[assignment]
                "mean": sum(s["mean"] for s in quality_stats_list) / len(quality_stats_list),
                "std": sum(s["std"] for s in quality_stats_list) / len(quality_stats_list)
            }

        return aggregated_stats

    def _aggregate_quality_metrics(self, results: list[dict[str, Any]]) -> dict[str, Any]:
        """Aggregate quality metrics from multiple chunks."""
        aggregated_metrics = {}

        total_records = sum(r.get("total_records", 0) for r in results)
        total_duplicates = sum(r.get("duplicate_records", 0) for r in results)

        aggregated_metrics["total_records"] = total_records
        aggregated_metrics["total_duplicate_records"] = total_duplicates

        # Aggregate missing values
        missing_values_combined = {}
        for result in results:
            if "missing_values" in result:
                for col, count in result["missing_values"].items():
                    missing_values_combined[col] = missing_values_combined.get(col, 0) + count

        aggregated_metrics["missing_values"] = missing_values_combined

        # Aggregate quality distribution
        quality_dist_combined = {
            "high_quality": 0,
            "medium_quality": 0,
            "low_quality": 0
        }

        for result in results:
            if "quality_distribution" in result:
                for category, count in result["quality_distribution"].items():
                    quality_dist_combined[category] += count

        aggregated_metrics["quality_distribution"] = quality_dist_combined

        return aggregated_metrics

    def _aggregate_filtering_results(self, results: list[dict[str, Any]]) -> dict[str, Any]:
        """Aggregate filtering results from multiple chunks."""
        aggregated_filtering = {}

        total_original = sum(r.get("original_count", 0) for r in results)
        total_filtered = sum(r.get("filtered_count", 0) for r in results)
        total_removed = sum(r.get("removed_count", 0) for r in results)

        aggregated_filtering["total_original_count"] = total_original
        aggregated_filtering["total_filtered_count"] = total_filtered
        aggregated_filtering["total_removed_count"] = total_removed
        aggregated_filtering["filter_efficiency"] = (total_removed / total_original) * 100 if total_original > 0 else 0

        # Get filter criteria from first result
        if results:
            aggregated_filtering["filter_criteria"] = results[0].get("filter_criteria", "Unknown")

        return aggregated_filtering

    def _aggregate_basic_results(self, results: list[dict[str, Any]]) -> dict[str, Any]:
        """Aggregate basic analysis results."""
        aggregated_basic = {}

        total_records = sum(r.get("record_count", 0) for r in results)
        total_memory_mb = sum(r.get("memory_usage_mb", 0) for r in results)

        aggregated_basic["total_record_count"] = total_records
        aggregated_basic["total_memory_usage_mb"] = total_memory_mb

        # Collect unique columns across all chunks
        all_columns = set()
        for result in results:
            if "columns" in result:
                all_columns.update(result["columns"])

        aggregated_basic["all_columns"] = sorted(list(all_columns))

        return aggregated_basic

    def _process_sequential(self, chunks: list[pd.DataFrame], analysis_type: str, **kwargs) -> dict[str, Any]:
        """Fallback sequential processing when Dask is not available."""
        results = []

        for i, chunk in enumerate(chunks):
            try:
                # Call the delayed function directly without Dask
                if analysis_type == "variant_stats":
                    result = self._compute_variant_statistics(chunk, **kwargs)
                elif analysis_type == "quality_metrics":
                    result = self._compute_quality_metrics(chunk, **kwargs)
                elif analysis_type == "annotation":
                    result = self._annotate_variants(chunk, **kwargs)
                elif analysis_type == "filtering":
                    result = self._filter_variants(chunk, **kwargs)
                else:
                    result = self._basic_analysis(chunk, **kwargs)

                result["chunk_size"] = len(chunk)
                results.append(result)
                logger.debug(f"Processed chunk {i+1}/{len(chunks)} sequentially")
            except Exception as e:
                logger.error(f"Error processing chunk {i}: {e}")
                results.append({"error": str(e)})

        return self._aggregate_results(results, analysis_type)

    def close(self) -> None:
        """Close Dask client and cleanup resources."""
        if self.client:
            self.client.close()
            self.client = None
            logger.info("Dask cluster closed")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
