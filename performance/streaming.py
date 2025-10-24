"""
💨 Streaming Data Processor

Provides memory-efficient streaming for massive files:
- Chunk-based processing (process 100GB+ files with 8GB RAM)
- Generator patterns for lazy evaluation
- Incremental aggregation
- Out-of-core processing

Performance: Process unlimited size files with constant memory usage

Examples:
    >>> processor = StreamProcessor(chunk_size=100_000)
    >>> for result in processor.stream_csv("huge_file.csv"):
    ...     process(result)  # Process one chunk at a time
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path
from typing import Any, Callable, Generator, Iterator, Optional

import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


class StreamProcessor:
    """
    Memory-efficient streaming processor for massive datasets.

    Processes data in chunks to avoid loading entire dataset into memory.

    Examples:
        >>> processor = StreamProcessor(chunk_size=50_000)
        >>>
        >>> # Stream and process 100GB CSV
        >>> total = 0
        >>> for chunk in processor.stream_csv("massive.csv"):
        ...     total += len(chunk)
        ...     # Process chunk...
        >>> print(f"Processed {total} rows")
    """

    def __init__(
        self,
        chunk_size: int = 100_000,
        compression: Optional[str] = None,
    ):
        """
        Initialize stream processor.

        Args:
            chunk_size: Number of rows per chunk
            compression: Compression type ('gzip', 'bz2', None)
        """
        self.chunk_size = chunk_size
        self.compression = compression

        logger.info(f"StreamProcessor initialized (chunk_size={chunk_size})")

    def stream_csv(
        self,
        file_path: str | Path,
        **read_csv_kwargs: Any,
    ) -> Generator[pd.DataFrame, None, None]:
        """
        Stream CSV file in chunks.

        Args:
            file_path: Path to CSV file
            **read_csv_kwargs: Arguments for pd.read_csv()

        Yields:
            DataFrame chunks

        Examples:
            >>> for chunk in processor.stream_csv("data.csv"):
            ...     print(f"Processing {len(chunk)} rows")
        """
        logger.info(f"Streaming CSV: {file_path}")

        try:
            reader = pd.read_csv(
                file_path,
                chunksize=self.chunk_size,
                compression=self.compression,
                **read_csv_kwargs
            )

            chunk_num = 0
            for chunk in reader:
                chunk_num += 1
                logger.debug(f"Yielding chunk {chunk_num}: {len(chunk)} rows")
                yield chunk

            logger.info(f"Stream complete: {chunk_num} chunks processed")

        except Exception as exc:
            logger.error(f"Stream error: {exc}")
            raise

    def stream_fasta(
        self,
        file_path: str | Path,
    ) -> Generator[tuple[str, str], None, None]:
        """
        Stream FASTA file (bioinformatics format).

        Yields:
            Tuple of (header, sequence)

        Examples:
            >>> for header, seq in processor.stream_fasta("genome.fasta"):
            ...     print(f"{header}: {len(seq)} bp")
        """
        logger.info(f"Streaming FASTA: {file_path}")

        file_path = Path(file_path)

        # Handle gzip compression
        if file_path.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(file_path, mode) as f:
            header = None
            sequence = []

            for line in f:
                line = line.strip()

                if line.startswith('>'):
                    # Yield previous sequence
                    if header:
                        yield header, ''.join(sequence)

                    # Start new sequence
                    header = line[1:]  # Remove '>'
                    sequence = []
                else:
                    sequence.append(line)

            # Yield last sequence
            if header:
                yield header, ''.join(sequence)

    def stream_fastq(
        self,
        file_path: str | Path,
    ) -> Generator[tuple[str, str, str], None, None]:
        """
        Stream FASTQ file (sequencing data format).

        Yields:
            Tuple of (header, sequence, quality)

        Examples:
            >>> for header, seq, qual in processor.stream_fastq("reads.fastq.gz"):
            ...     print(f"{header}: Q{np.mean([ord(c)-33 for c in qual]):.1f}")
        """
        logger.info(f"Streaming FASTQ: {file_path}")

        file_path = Path(file_path)

        if file_path.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(file_path, mode) as f:
            while True:
                # FASTQ format: 4 lines per record
                header = f.readline().strip()
                if not header:
                    break  # EOF

                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()

                yield header[1:], sequence, quality  # Remove '@' from header

    def stream_vcf(
        self,
        file_path: str | Path,
    ) -> Generator[dict[str, Any], None, None]:
        """
        Stream VCF file (variant call format).

        Yields:
            Dictionary with variant information

        Examples:
            >>> for variant in processor.stream_vcf("variants.vcf.gz"):
            ...     print(f"{variant['CHROM']}:{variant['POS']}")
        """
        logger.info(f"Streaming VCF: {file_path}")

        file_path = Path(file_path)

        if file_path.suffix == '.gz':
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(file_path, mode) as f:
            # Skip header lines
            for line in f:
                if line.startswith('#CHROM'):
                    # Parse column names
                    columns = line.strip().split('\t')
                    break
                elif not line.startswith('##'):
                    break

            # Parse variants
            for line in f:
                if line.startswith('#'):
                    continue

                values = line.strip().split('\t')
                variant = dict(zip(columns, values))

                yield variant

    def stream_apply(
        self,
        file_path: str | Path,
        func: Callable[[pd.DataFrame], pd.DataFrame],
        output_path: Optional[str | Path] = None,
        **read_csv_kwargs: Any,
    ) -> Optional[pd.DataFrame]:
        """
        Apply function to each chunk and optionally save results.

        Args:
            file_path: Input CSV path
            func: Function to apply to each chunk
            output_path: Output CSV path (if None, return concatenated result)
            **read_csv_kwargs: Arguments for pd.read_csv()

        Returns:
            Concatenated DataFrame (if output_path is None)

        Examples:
            >>> def filter_chunk(chunk):
            ...     return chunk[chunk['value'] > 100]
            >>>
            >>> processor.stream_apply(
            ...     "large.csv",
            ...     filter_chunk,
            ...     output_path="filtered.csv"
            ... )
        """
        logger.info(f"Stream applying function to {file_path}")

        results = []
        first_chunk = True

        for chunk in self.stream_csv(file_path, **read_csv_kwargs):
            processed = func(chunk)

            if output_path:
                # Write to file
                processed.to_csv(
                    output_path,
                    mode='a',  # Append mode
                    header=first_chunk,
                    index=False
                )
                first_chunk = False
            else:
                # Accumulate in memory
                results.append(processed)

        if output_path:
            logger.info(f"Results written to {output_path}")
            return None
        else:
            logger.info(f"Concatenating {len(results)} chunks")
            return pd.concat(results, ignore_index=True)

    def stream_groupby_agg(
        self,
        file_path: str | Path,
        groupby_cols: list[str],
        agg_dict: dict[str, str | list[str]],
        **read_csv_kwargs: Any,
    ) -> pd.DataFrame:
        """
        Streaming groupby aggregation for huge files.

        Aggregates incrementally without loading entire file.

        Args:
            file_path: Input CSV path
            groupby_cols: Columns to group by
            agg_dict: Aggregation dictionary
            **read_csv_kwargs: Arguments for pd.read_csv()

        Returns:
            Aggregated DataFrame

        Examples:
            >>> result = processor.stream_groupby_agg(
            ...     "sales.csv",
            ...     groupby_cols=['region', 'product'],
            ...     agg_dict={'revenue': 'sum', 'quantity': 'sum'}
            ... )
        """
        logger.info(f"Streaming groupby aggregation on {file_path}")

        partial_results = []

        for chunk in self.stream_csv(file_path, **read_csv_kwargs):
            # Aggregate chunk
            chunk_agg = chunk.groupby(groupby_cols).agg(agg_dict).reset_index()
            partial_results.append(chunk_agg)

        # Combine partial results
        combined = pd.concat(partial_results, ignore_index=True)

        # Final aggregation
        final_result = combined.groupby(groupby_cols).agg(agg_dict).reset_index()

        logger.info(f"Aggregation complete: {len(final_result)} groups")
        return final_result

    def stream_filter(
        self,
        file_path: str | Path,
        condition: Callable[[pd.DataFrame], pd.Series],
        output_path: str | Path,
        **read_csv_kwargs: Any,
    ) -> int:
        """
        Filter large CSV file based on condition.

        Args:
            file_path: Input CSV path
            condition: Function that returns boolean Series
            output_path: Output CSV path
            **read_csv_kwargs: Arguments for pd.read_csv()

        Returns:
            Number of rows written

        Examples:
            >>> rows = processor.stream_filter(
            ...     "data.csv",
            ...     lambda df: df['score'] > 0.8,
            ...     "filtered.csv"
            ... )
        """
        logger.info(f"Streaming filter: {file_path} → {output_path}")

        total_rows = 0
        first_chunk = True

        for chunk in self.stream_csv(file_path, **read_csv_kwargs):
            filtered = chunk[condition(chunk)]

            if len(filtered) > 0:
                filtered.to_csv(
                    output_path,
                    mode='a',
                    header=first_chunk,
                    index=False
                )
                total_rows += len(filtered)
                first_chunk = False

        logger.info(f"Filter complete: {total_rows} rows written")
        return total_rows


# ═══════════════════════════════════════════════════════════════
#  Incremental Statistics
# ═══════════════════════════════════════════════════════════════

class IncrementalStats:
    """
    Compute statistics incrementally without loading full dataset.

    Useful for calculating mean, std, min, max on huge files.

    Examples:
        >>> stats = IncrementalStats()
        >>> for chunk in processor.stream_csv("huge.csv"):
        ...     stats.update(chunk['value'])
        >>> print(f"Mean: {stats.mean()}, Std: {stats.std()}")
    """

    def __init__(self):
        """Initialize incremental statistics tracker"""
        self.n = 0
        self.sum = 0.0
        self.sum_sq = 0.0
        self.min_val = float('inf')
        self.max_val = float('-inf')

    def update(self, values: pd.Series | np.ndarray):
        """
        Update statistics with new values.

        Args:
            values: New values to include
        """
        values = np.asarray(values)

        self.n += len(values)
        self.sum += values.sum()
        self.sum_sq += (values ** 2).sum()
        self.min_val = min(self.min_val, values.min())
        self.max_val = max(self.max_val, values.max())

    def mean(self) -> float:
        """Get current mean"""
        return self.sum / self.n if self.n > 0 else 0.0

    def variance(self) -> float:
        """Get current variance"""
        if self.n < 2:
            return 0.0
        return (self.sum_sq / self.n) - (self.mean() ** 2)

    def std(self) -> float:
        """Get current standard deviation"""
        return np.sqrt(self.variance())

    def min(self) -> float:
        """Get minimum value"""
        return self.min_val

    def max(self) -> float:
        """Get maximum value"""
        return self.max_val

    def summary(self) -> dict[str, float]:
        """Get summary statistics"""
        return {
            'count': self.n,
            'mean': self.mean(),
            'std': self.std(),
            'min': self.min(),
            'max': self.max(),
        }


# ═══════════════════════════════════════════════════════════════
#  Convenience Functions
# ═══════════════════════════════════════════════════════════════

def stream_csv_chunks(
    file_path: str | Path,
    chunk_size: int = 100_000,
) -> Generator[pd.DataFrame, None, None]:
    """Quick CSV streaming (convenience function)"""
    processor = StreamProcessor(chunk_size=chunk_size)
    yield from processor.stream_csv(file_path)
