"""
GenoScope main entry point with improved error handling and structure.
"""

import argparse
import logging
import sys
from pathlib import Path

import pandas as pd
from sklearn.preprocessing import LabelEncoder

from genoscope.core.logging_config import setup_logging
from genoscope.data_analysis.analysis_core import extract_pca
from genoscope.data_analysis.data_cleaning import handle_missing_values
from genoscope.data_analysis.data_cleaning import remove_duplicates
from genoscope.data_analysis.data_filtering import filter_by_custom_function
from genoscope.data_analysis.data_filtering import filter_by_multiple_conditions
from genoscope.data_analysis.data_filtering import filter_by_percentile
from genoscope.data_analysis.data_filtering import filter_outliers
from genoscope.data_analysis.data_ingestion import load_data
from genoscope.data_analysis.visualization import plot_pca

from typing import Dict, List, Tuple, Set, Optional, Union, Any
# Setup logging
logger = logging.getLogger(__name__)
setup_logging()


class GenoScopeProcessor:
    """Main processor class for GenoScope data analysis pipeline."""

    def __init__(self, config: Optional[dict] = None):
        self.config = config or {}
        self.data: pd.DataFrame | None = None

        # Добавляем поддержку параллелизации
        self._parallel_enabled = False
        self._dask_processor = None

    def set_parallel_config(self,
                          enable: bool = True,
                          n_workers: int = 4,
                          memory_limit: str = "2GB") -> None:
        """Настройка параллельной обработки.
        
        Args:
            enable: Включить параллельную обработку
            n_workers: Количество воркеров
            memory_limit: Лимит памяти на воркер
        """
        self._parallel_enabled = enable

        if enable:
            try:
                from genoscope.parallel import DaskGenomicProcessor
                self._dask_processor = DaskGenomicProcessor(
                    n_workers=n_workers,
                    memory_limit=memory_limit
                )
                logger.info(f"Parallel processing enabled with {n_workers} workers")
            except ImportError:
                logger.warning("Dask not available, parallel processing disabled")
                self._parallel_enabled = False

    def load_data_enhanced(self,
                          file_path: str,
                          file_type: str,
                          force_parallel: bool = False) -> bool:
        """Улучшенная загрузка с поддержкой параллелизации.
        
        Args:
            file_path: Путь к файлу
            file_type: Тип файла
            force_parallel: Принудительно использовать параллельную обработку
            
        Returns:
            bool: True если успешно
        """

        if self._parallel_enabled and (force_parallel or self._should_use_parallel(file_path)):
            return self._load_parallel(file_path, file_type)
        return self.load_file(file_path, file_type)

    def _should_use_parallel(self, file_path: str) -> bool:
        """Определяет, нужна ли параллельная обработка."""
        try:
            file_size_mb = Path(file_path).stat().st_size / (1024 * 1024)
            return file_size_mb > 100  # Файлы больше 100MB обрабатываем параллельно
        except Exception:
            return False

    def _load_parallel(self, file_path: str, file_type: str) -> bool:
        """Параллельная загрузка больших файлов."""
        try:
            results = self._dask_processor.process_large_file_parallel(
                file_path=file_path,
                file_type=file_type,
                analysis_type="basic"
            )

            if "error" not in results:
                perf_summary = results.get("performance_summary", {})
                total_records = perf_summary.get("total_records_processed", 0)
                processing_speed = perf_summary.get("average_processing_speed", 0)

                logger.info(f"Parallel processing completed: {total_records:,} records at {processing_speed:.0f} records/sec")

                # Сохраняем результаты для дальнейшей обработки
                self._parallel_results = results
                return True
            logger.error(f"Parallel processing failed: {results['error']}")
            return False

        except Exception as e:
            logger.error(f"Parallel processing failed: {e}")
            logger.info("Falling back to sequential processing")
            return self.load_file(file_path, file_type)

    def run_parallel_analysis(self, analysis_type: str = "comprehensive") -> dict:
        """Запустить параллельный анализ.
        
        Args:
            analysis_type: Тип анализа (basic, variant_stats, quality_metrics, etc.)
            
        Returns:
            Результаты анализа
        """
        if not self._parallel_enabled or not self._dask_processor:
            raise RuntimeError("Parallel processing not enabled")

        logger.info(f"Running parallel {analysis_type} analysis")

        # Здесь могут быть различные типы анализа
        # Пока возвращаем базовую информацию
        return {
            "analysis_type": analysis_type,
            "status": "completed",
            "parallel_enabled": True
        }

    def get_performance_report(self) -> dict:
        """Получить отчет о производительности."""
        if not self._dask_processor:
            return {"error": "Parallel processing not enabled"}

        if hasattr(self, "_parallel_results"):
            return self._parallel_results.get("performance", {})
        return {"status": "No parallel operations performed yet"}

    def export_performance_metrics(self, output_path: str, format: str = "html"):
        """Экспортировать метрики производительности.
        
        Args:
            output_path: Путь к выходному файлу
            format: Формат отчета (html, json)
            
        Returns:
            Путь к созданному файлу
        """
        if not self._dask_processor:
            raise RuntimeError("Parallel processing not enabled")

        return self._dask_processor.performance_monitor.export_report(
            Path(output_path), format
        )

    def load_file(self, file_path: str, file_type: str) -> bool:
        """
        Load data from file with error handling.

        Args:
            file_path: Path to the data file
            file_type: Type of the file (csv, vcf, etc.)

        Returns:
            bool: True if successful, False otherwise
        """
        if not Path(file_path).exists():
            logger.error(f"File not found: {file_path}")
            return False

        try:
            self.data = load_data(file_path, file_type)
            if self.data is None:
                logger.error("Failed to load data - returned None")
                return False

            logger.info(f"Successfully loaded {len(self.data)} rows from {file_path}")
            return True

        except Exception as e:
            logger.error(f"Error loading data from {file_path}: {e}")
            return False

    def clean_data(self) -> bool:
        """
        Clean the loaded data.

        Returns:
            bool: True if successful, False otherwise
        """
        if self.data is None:
            logger.error("No data loaded for cleaning")
            return False

        try:
            # Remove duplicates
            original_size = len(self.data)
            self.data = remove_duplicates(self.data)
            logger.info(f"Removed {original_size - len(self.data)} duplicate rows")

            # Handle missing values
            missing_count = self.data.isna().sum().sum()
            if missing_count > 0:
                self.data = handle_missing_values(self.data, method="mean")
                logger.info(f"Handled {missing_count} missing values")

            # Encode categorical data
            categorical_cols = self.data.select_dtypes(
                include=["object", "category"]
            ).columns
            for column in categorical_cols:
                try:
                    le = LabelEncoder()
                    self.data[column] = le.fit_transform(self.data[column])
                    logger.info(f"Encoded categorical column: {column}")
                except Exception as e:
                    logger.warning(f"Failed to encode column {column}: {e}")

            return True

        except Exception as e:
            logger.error(f"Error during data cleaning: {e}")
            return False

    def apply_filters(self) -> bool:
        """
        Apply various filters to the data.

        Returns:
            bool: True if successful, False otherwise
        """
        if self.data is None:
            logger.error("No data loaded for filtering")
            return False

        try:
            original_size = len(self.data)

            # Apply multiple condition filters if numeric columns exist
            numeric_cols = self.data.select_dtypes(include=["float64", "int"]).columns
            if len(numeric_cols) >= 2:
                first_col, second_col = numeric_cols[0], numeric_cols[1]

                # Dynamic conditions based on data quartiles
                q1_first = self.data[first_col].quantile(0.25)
                q3_second = self.data[second_col].quantile(0.75)

                conditions = [
                    f"{first_col} > {q1_first}",
                    f"{second_col} < {q3_second}",
                ]
                self.data = filter_by_multiple_conditions(self.data, conditions)
                logger.info(f"Applied condition filters: {conditions}")

                # Apply custom function filter
                self.data = filter_by_custom_function(
                    self.data,
                    lambda row: row[first_col] > self.data[first_col].median(),
                )
                logger.info("Applied custom function filter")

                # Apply percentile filter
                self.data = filter_by_percentile(
                    self.data, first_col, lower=10, upper=90
                )
                logger.info(f"Applied percentile filter to {first_col}")

                # Apply outlier filter
                self.data = filter_outliers(self.data, first_col, method="iqr")
                logger.info(f"Applied outlier filter to {first_col}")

            filtered_size = len(self.data)
            logger.info(
                f"Filtering reduced data from {original_size} to {filtered_size} rows"
            )
            return True

        except Exception as e:
            logger.error(f"Error during filtering: {e}")
            return False

    def perform_pca(self, n_components: int = 2) -> pd.DataFrame | None:
        """
        Perform PCA analysis on numeric data.

        Args:
            n_components: Number of principal components

        Returns:
            DataFrame with principal components or None if failed
        """
        if self.data is None:
            logger.error("No data loaded for PCA")
            return None

        try:
            # Select only numeric columns
            numeric_data = self.data.select_dtypes(include=["float64", "int"])
            if numeric_data.empty:
                logger.error("No numeric columns found for PCA")
                return None

            if len(numeric_data.columns) < n_components:
                logger.warning(
                    f"Only {len(numeric_data.columns)} numeric columns, reducing components"
                )
                n_components = len(numeric_data.columns)

            pca_result = extract_pca(numeric_data, n_components=n_components)
            logger.info(f"PCA completed with {n_components} components")

            # Try to plot if visualization is available
            try:
                target_column = None  # Could be enhanced to detect target column
                plot_pca(pca_result, target_column)
                logger.info("PCA plot generated")
            except Exception as plot_error:
                logger.warning(f"Could not generate PCA plot: {plot_error}")

            return pca_result

        except Exception as e:
            logger.error(f"Error during PCA analysis: {e}")
            return None

    def run_pipeline(
        self, file_path: str, file_type: str, output_path: Optional[str] = None
    ) -> bool:
        """
        Run the complete analysis pipeline.

        Args:
            file_path: Path to input file
            file_type: Type of input file
            output_path: Optional output path for results

        Returns:
            bool: True if pipeline completed successfully
        """
        logger.info(f"Starting GenoScope pipeline for {file_path}")

        # Step 1: Load data
        if not self.load_file(file_path, file_type):
            return False

        # Step 2: Clean data
        if not self.clean_data():
            return False

        # Step 3: Apply filters
        if not self.apply_filters():
            logger.warning("Filtering step failed, continuing with unfiltered data")

        # Step 4: PCA analysis
        pca_result = self.perform_pca()
        if pca_result is not None:
            logger.info(f"PCA result shape: {pca_result.shape}")

        # Step 5: Save results if output path specified
        if output_path and self.data is not None:
            try:
                output_file = Path(output_path)
                output_file.parent.mkdir(parents=True, exist_ok=True)
                self.data.to_csv(output_file, index=False)
                logger.info(f"Results saved to {output_path}")
            except Exception as e:
                logger.error(f"Failed to save results: {e}")

        logger.info("Pipeline completed successfully")
        return True


# GUI functionality removed


def main():
    """Main entry point with CLI argument parsing."""
    parser = argparse.ArgumentParser(
        description="GenoScope - Genomic data analysis toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --input data.csv --type csv       # Analyze CSV file
  %(prog)s --input data.vcf --type vcf --output result.csv  # VCF to CSV
  %(prog)s --input large.csv --type csv --parallel --workers 8  # Parallel processing
        """,
    )

    parser.add_argument("--input", "-i", type=str, help="Input file path")
    parser.add_argument(
        "--type", "-t", type=str, help="File type (csv, vcf, bam, etc.)"
    )
    parser.add_argument("--output", "-o", type=str, help="Output file path")
# GUI option removed
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")

    # Parallel processing options
    parser.add_argument("--parallel", action="store_true", help="Enable parallel processing")
    parser.add_argument("--workers", type=int, default=4, help="Number of parallel workers (default: 4)")
    parser.add_argument("--memory-limit", type=str, default="2GB", help="Memory limit per worker (default: 2GB)")
    parser.add_argument("--analysis-type", type=str, default="basic",
                       choices=["basic", "variant_stats", "quality_metrics", "annotation", "filtering"],
                       help="Type of analysis to perform")
    parser.add_argument("--performance-report", type=str,
                       help="Export performance report to file (HTML format)")

    args = parser.parse_args()

    # Set log level based on verbose flag
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        if args.input and args.type:
            processor = GenoScopeProcessor()

            # Configure parallel processing if requested
            if args.parallel:
                processor.set_parallel_config(
                    enable=True,
                    n_workers=args.workers,
                    memory_limit=args.memory_limit
                )
                logger.info(f"Parallel processing enabled: {args.workers} workers, {args.memory_limit} memory limit")

            # Use enhanced loading if parallel processing is enabled
            if args.parallel:
                success = processor.load_data_enhanced(args.input, args.type, force_parallel=True)
            else:
                success = processor.run_pipeline(args.input, args.type, args.output)

            # Export performance report if requested
            if args.performance_report and args.parallel:
                try:
                    report_path = processor.export_performance_metrics(args.performance_report, "html")
                    logger.info(f"Performance report exported to: {report_path}")
                except Exception as e:
                    logger.warning(f"Could not export performance report: {e}")

            # Run specific analysis if requested and parallel processing is enabled
            if args.parallel and success:
                try:
                    analysis_results = processor.run_parallel_analysis(args.analysis_type)
                    logger.info(f"Parallel {args.analysis_type} analysis completed: {analysis_results}")
                except Exception as e:
                    logger.warning(f"Parallel analysis failed: {e}")

            sys.exit(0 if success else 1)
        else:
            # Show help if no arguments provided
            logger.info("No arguments provided. Use --help for usage information.")
            parser.print_help()
            sys.exit(1)

    except KeyboardInterrupt:
        logger.info("Operation cancelled by user")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
