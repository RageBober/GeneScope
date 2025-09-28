"""Intelligent chunk managers for parallel processing of genomic data.

This module provides smart chunking strategies for different types of genomic files,
optimizing for both performance and biological relevance.
"""

from pathlib import Path

import pandas as pd

from ..core.logging_config import get_logger

logger = get_logger(__name__)


class BaseChunkManager:
    """Base class for file chunking strategies."""

    def __init__(self):
        self.logger = logger

    def create_chunks(self, file_path: Path, chunk_size_mb: int = 50) -> list[pd.DataFrame]:
        """Create chunks from file.
        
        Args:
            file_path: Path to input file
            chunk_size_mb: Target size for each chunk in MB
            
        Returns:
            List of DataFrame chunks
        """
        raise NotImplementedError("Subclasses must implement create_chunks")

    def estimate_chunk_size(self, file_path: Path, target_size_mb: int) -> int:
        """Estimate optimal chunk size in rows.
        
        Args:
            file_path: Path to input file
            target_size_mb: Target chunk size in MB
            
        Returns:
            Estimated number of rows per chunk
        """
        try:
            # Read sample to estimate memory usage
            sample_size = min(1000, self._count_lines(file_path) - 1)
            if sample_size <= 0:
                return 1000  # fallback

            sample_df = self._read_sample(file_path, sample_size)
            bytes_per_row = sample_df.memory_usage(deep=True).sum() / len(sample_df)

            target_bytes = target_size_mb * 1024 * 1024
            rows_per_chunk = int(target_bytes / bytes_per_row)

            # Ensure reasonable bounds
            return max(100, min(rows_per_chunk, 100000))

        except Exception as e:
            self.logger.warning(f"Could not estimate chunk size: {e}")
            return 10000  # reasonable fallback

    def _count_lines(self, file_path: Path) -> int:
        """Count lines in file efficiently."""
        try:
            with open(file_path, "rb") as f:
                lines = sum(1 for _ in f)
            return lines
        except Exception:
            return 1000  # fallback

    def _read_sample(self, file_path: Path, nrows: int) -> pd.DataFrame:
        """Read sample from file for estimation."""
        return pd.read_csv(file_path, nrows=nrows)

    def _estimate_size_mb(self, df: pd.DataFrame) -> float:
        """Estimate DataFrame size in MB."""
        return float(df.memory_usage(deep=True).sum() / (1024 * 1024))


class CSVChunkManager(BaseChunkManager):
    """Smart chunk manager for CSV files with automatic delimiter detection."""

    def __init__(self):
        super().__init__()
        self._delimiter = None

    def create_chunks(self, file_path: Path, chunk_size_mb: int = 50) -> list[pd.DataFrame]:
        """Create optimally-sized chunks from CSV file.
        
        Args:
            file_path: Path to CSV file
            chunk_size_mb: Target chunk size in MB
            
        Returns:
            List of DataFrame chunks
        """
        try:
            # Detect delimiter and estimate chunk size
            self._detect_separator(file_path)
            chunk_size_rows = self.estimate_chunk_size(file_path, chunk_size_mb)

            self.logger.info(f"Processing CSV with delimiter='{self._delimiter}', chunk_size={chunk_size_rows} rows")

            chunks = []
            chunk_reader = pd.read_csv(
                file_path,
                sep=self._delimiter,
                chunksize=chunk_size_rows,
                low_memory=False
            )

            for i, chunk_df in enumerate(chunk_reader):
                if not chunk_df.empty:
                    chunks.append(chunk_df)
                    self.logger.debug(f"Created chunk {i+1}: {len(chunk_df)} rows")

            self.logger.info(f"Created {len(chunks)} chunks from {file_path.name}")
            return chunks

        except Exception as e:
            self.logger.error(f"Failed to chunk CSV file {file_path}: {e}")
            raise

    def _detect_separator(self, file_path: Path) -> str:
        """Auto-detect CSV separator/delimiter."""
        if self._delimiter:
            return str(self._delimiter)

        try:
            # Read first few lines to detect separator
            with open(file_path, encoding="utf-8") as f:
                first_lines = [f.readline() for _ in range(5)]

            # Try common separators
            separators = [",", "\t", ";", "|"]
            separator_scores = {}

            for sep in separators:
                try:
                    df_test = pd.read_csv(file_path, sep=sep, nrows=10)
                    # Score based on number of columns and consistency
                    score = len(df_test.columns) * (1.0 if len(df_test) > 0 else 0.5)
                    separator_scores[sep] = score
                except Exception:
                    separator_scores[sep] = 0

            # Choose best separator
            best_sep = str(max(separator_scores.items(), key=lambda x: x[1])[0])
            self._delimiter = best_sep

            self.logger.debug(f"Detected CSV separator: '{self._delimiter}'")
            return self._delimiter

        except Exception as e:
            self.logger.warning(f"Could not detect separator, using comma: {e}")
            self._delimiter = ","
            return self._delimiter

    def _read_sample(self, file_path: Path, nrows: int) -> pd.DataFrame:
        """Read sample from CSV with detected separator."""
        return pd.read_csv(file_path, sep=self._delimiter, nrows=nrows)


class VCFChunkManager(BaseChunkManager):
    """Biologically-aware chunk manager for VCF files."""

    def create_chunks(self, file_path: Path, chunk_size_mb: int = 50) -> list[pd.DataFrame]:
        """Create biologically-meaningful chunks from VCF file.
        
        VCF chunking strategy:
        1. Split by chromosome (biologically meaningful)
        2. If chromosome too large, split by genomic position
        3. Maintain header information consistency
        
        Args:
            file_path: Path to VCF file
            chunk_size_mb: Target chunk size in MB
            
        Returns:
            List of DataFrame chunks
        """
        try:
            # Read VCF file (skip comment lines)
            full_df = self._read_vcf(file_path)

            if full_df.empty:
                self.logger.warning(f"VCF file {file_path} is empty")
                return []

            # Strategy 1: Chunk by chromosome
            if "CHROM" in full_df.columns:
                chunks = self._chunk_by_chromosome(full_df, chunk_size_mb)
            else:
                # Fallback: simple size-based chunking
                self.logger.warning("No CHROM column found, using simple chunking")
                chunks = self._simple_chunk(full_df, chunk_size_mb)

            self.logger.info(f"Created {len(chunks)} biologically-aware chunks from VCF")
            return chunks

        except Exception as e:
            self.logger.error(f"Failed to chunk VCF file {file_path}: {e}")
            raise

    def _read_vcf(self, file_path: Path) -> pd.DataFrame:
        """Read VCF file, handling header and comment lines."""
        try:
            # Find header line (starts with #CHROM)
            header_line = None
            with open(file_path) as f:
                for i, line in enumerate(f):
                    if line.startswith("#CHROM"):
                        header_line = i
                        break

            if header_line is not None:
                # Read VCF with proper header
                df = pd.read_csv(
                    file_path,
                    sep="\t",
                    skiprows=header_line,
                    comment="#",
                    low_memory=False
                )
                # Clean column names (remove leading #)
                df.columns = [col.lstrip("#") for col in df.columns]
            else:
                # Try reading as regular tab-separated file
                df = pd.read_csv(file_path, sep="\t", comment="#", low_memory=False)

            self.logger.debug(f"Read VCF with {len(df)} variants and {len(df.columns)} columns")
            return df

        except Exception as e:
            self.logger.error(f"Failed to read VCF file: {e}")
            raise

    def _chunk_by_chromosome(self, df: pd.DataFrame, chunk_size_mb: int) -> list[pd.DataFrame]:
        """Chunk VCF data by chromosome."""
        chunks = []

        for chrom in sorted(df["CHROM"].unique()):
            chrom_data = df[df["CHROM"] == chrom].copy()

            # Check if this chromosome fits in target size
            chrom_size_mb = self._estimate_size_mb(chrom_data)

            if chrom_size_mb <= chunk_size_mb:
                # Chromosome fits in one chunk
                chunks.append(chrom_data)
                self.logger.debug(f"Chromosome {chrom}: {len(chrom_data)} variants in one chunk")
            else:
                # Need to split chromosome by position
                sub_chunks = self._split_by_position(chrom_data, chunk_size_mb)
                chunks.extend(sub_chunks)
                self.logger.debug(f"Chromosome {chrom}: split into {len(sub_chunks)} chunks")

        return chunks

    def _split_by_position(self, chrom_df: pd.DataFrame, chunk_size_mb: int) -> list[pd.DataFrame]:
        """Split chromosome data by genomic position."""
        if "POS" not in chrom_df.columns:
            # Fallback to simple chunking
            return self._simple_chunk(chrom_df, chunk_size_mb)

        # Sort by position
        chrom_df = chrom_df.sort_values("POS")

        chunks = []
        chunk_size_rows = max(1000, len(chrom_df) // 10)  # At least 1000 variants per chunk

        for start_idx in range(0, len(chrom_df), chunk_size_rows):
            end_idx = min(start_idx + chunk_size_rows, len(chrom_df))
            chunk = chrom_df.iloc[start_idx:end_idx].copy()

            if not chunk.empty:
                chunks.append(chunk)

        return chunks

    def _simple_chunk(self, df: pd.DataFrame, chunk_size_mb: int) -> list[pd.DataFrame]:
        """Simple row-based chunking fallback."""
        # Use a fallback estimation based on average row size
        avg_row_size = max(1, int(self._estimate_size_mb(df) * 1024 * 1024 / len(df))) if len(df) > 0 else 1000
        estimated_rows = max(1, int(chunk_size_mb * 1024 * 1024 / avg_row_size))
        chunks = []

        for start_idx in range(0, len(df), estimated_rows):
            end_idx = min(start_idx + estimated_rows, len(df))
            chunk = df.iloc[start_idx:end_idx].copy()

            if not chunk.empty:
                chunks.append(chunk)

        return chunks


# Factory function for easy chunk manager creation
def get_chunk_manager(file_type: str) -> BaseChunkManager:
    """Get appropriate chunk manager for file type.
    
    Args:
        file_type: Type of file ('csv', 'vcf', etc.)
        
    Returns:
        Appropriate chunk manager instance
    """
    file_type = file_type.lower()

    if file_type in ["csv", "tsv"]:
        return CSVChunkManager()
    if file_type == "vcf":
        return VCFChunkManager()
    # Default to CSV manager for unknown types
    logger.warning(f"Unknown file type '{file_type}', using CSV chunk manager")
    return CSVChunkManager()
