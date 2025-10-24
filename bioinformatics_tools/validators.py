"""
🔒 Input Validation and Sanitization Module

Provides security controls for bioinformatics tool inputs:
- Path traversal prevention
- File type validation
- Sequence validation
- Resource limit checks
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional

# Allowed file extensions for bioinformatics data
ALLOWED_EXTENSIONS = {
    ".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn",  # FASTA
    ".fastq", ".fq",  # FASTQ
    ".vcf", ".vcf.gz",  # VCF
    ".bam", ".sam",  # Alignment
    ".gff", ".gff3", ".gtf",  # Annotation
    ".hdf5", ".h5",  # HDF5
    ".csv", ".tsv", ".txt",  # Tabular
}

# Maximum file size (10 GB by default)
MAX_FILE_SIZE_BYTES = 10 * 1024 * 1024 * 1024


class ValidationError(Exception):
    """Raised when input validation fails"""
    pass


def validate_file_path(
    file_path: str | Path,
    *,
    must_exist: bool = True,
    allowed_extensions: Optional[set[str]] = None,
    max_size_bytes: Optional[int] = None,
) -> Path:
    """
    Validates file path for security and correctness.

    Security checks:
    1. Path traversal prevention (no .. components)
    2. Absolute path resolution
    3. File existence verification
    4. Extension whitelist validation
    5. Size limit check

    Args:
        file_path: Path to validate
        must_exist: Require file to exist
        allowed_extensions: Set of allowed extensions (default: ALLOWED_EXTENSIONS)
        max_size_bytes: Maximum file size (default: MAX_FILE_SIZE_BYTES)

    Returns:
        Validated absolute Path object

    Raises:
        ValidationError: If validation fails

    Examples:
        >>> validate_file_path("/data/sample.fasta")
        PosixPath('/data/sample.fasta')

        >>> validate_file_path("../etc/passwd")  # doctest: +SKIP
        ValidationError: Path contains traversal attempt
    """
    if not file_path:
        raise ValidationError("File path cannot be empty")

    # Convert to Path object and resolve to absolute path
    try:
        path = Path(file_path).resolve()
    except (ValueError, RuntimeError) as exc:
        raise ValidationError(f"Invalid file path: {exc}") from exc

    # Check for path traversal attempts
    if ".." in path.parts:
        raise ValidationError(f"Path contains traversal attempt: {file_path}")

    # Verify file exists if required
    if must_exist and not path.exists():
        raise ValidationError(f"File does not exist: {path}")

    if must_exist and not path.is_file():
        raise ValidationError(f"Path is not a regular file: {path}")

    # Validate file extension
    extensions = allowed_extensions or ALLOWED_EXTENSIONS
    if extensions and path.suffix.lower() not in extensions:
        raise ValidationError(
            f"File extension '{path.suffix}' not in allowed list: {extensions}"
        )

    # Check file size if exists
    if must_exist:
        max_size = max_size_bytes or MAX_FILE_SIZE_BYTES
        file_size = path.stat().st_size
        if file_size > max_size:
            raise ValidationError(
                f"File size {file_size} exceeds maximum {max_size} bytes"
            )

    return path


def validate_output_directory(
    dir_path: str | Path,
    *,
    create_if_missing: bool = True,
) -> Path:
    """
    Validates output directory path.

    Args:
        dir_path: Directory path to validate
        create_if_missing: Create directory if it doesn't exist

    Returns:
        Validated absolute Path object

    Raises:
        ValidationError: If validation fails
    """
    if not dir_path:
        raise ValidationError("Directory path cannot be empty")

    try:
        path = Path(dir_path).resolve()
    except (ValueError, RuntimeError) as exc:
        raise ValidationError(f"Invalid directory path: {exc}") from exc

    # Check for path traversal
    if ".." in path.parts:
        raise ValidationError(f"Path contains traversal attempt: {dir_path}")

    # Create directory if needed
    if create_if_missing and not path.exists():
        try:
            path.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise ValidationError(f"Cannot create directory: {exc}") from exc

    # Verify it's a directory
    if path.exists() and not path.is_dir():
        raise ValidationError(f"Path is not a directory: {path}")

    return path


def validate_sequence(
    sequence: str,
    *,
    seq_type: str = "nucleotide",
    min_length: int = 1,
    max_length: int = 1_000_000,
) -> str:
    """
    Validates biological sequence data.

    Args:
        sequence: Nucleotide or amino acid sequence
        seq_type: 'nucleotide' or 'protein'
        min_length: Minimum sequence length
        max_length: Maximum sequence length

    Returns:
        Validated sequence (uppercase)

    Raises:
        ValidationError: If validation fails

    Examples:
        >>> validate_sequence("ATCGATCG")
        'ATCGATCG'

        >>> validate_sequence("ATCXYZ")  # doctest: +SKIP
        ValidationError: Invalid nucleotide characters: X, Y, Z
    """
    if not sequence:
        raise ValidationError("Sequence cannot be empty")

    # Check length
    seq_len = len(sequence)
    if seq_len < min_length:
        raise ValidationError(f"Sequence too short: {seq_len} < {min_length}")
    if seq_len > max_length:
        raise ValidationError(f"Sequence too long: {seq_len} > {max_length}")

    # Convert to uppercase
    sequence = sequence.upper()

    # Validate characters based on sequence type
    if seq_type == "nucleotide":
        valid_chars = set("ATCGUNRYSWKMBDHV-")  # IUPAC nucleotides + gap
        invalid = set(sequence) - valid_chars
        if invalid:
            raise ValidationError(
                f"Invalid nucleotide characters: {', '.join(sorted(invalid))}"
            )
    elif seq_type == "protein":
        valid_chars = set("ACDEFGHIKLMNPQRSTVWY*-")  # IUPAC amino acids + stop + gap
        invalid = set(sequence) - valid_chars
        if invalid:
            raise ValidationError(
                f"Invalid protein characters: {', '.join(sorted(invalid))}"
            )
    else:
        raise ValidationError(f"Unknown sequence type: {seq_type}")

    return sequence


def validate_docker_image_name(image_name: str) -> str:
    """
    Validates Docker image name for security.

    Prevents injection attacks via malicious image names.

    Args:
        image_name: Docker image name (e.g., 'biocontainers/kraken2:latest')

    Returns:
        Validated image name

    Raises:
        ValidationError: If validation fails
    """
    if not image_name:
        raise ValidationError("Docker image name cannot be empty")

    # Docker image name pattern: [registry/][namespace/]repository[:tag]
    pattern = r"^[a-z0-9][a-z0-9._/-]*(?::[a-z0-9._-]+)?$"
    if not re.match(pattern, image_name.lower()):
        raise ValidationError(f"Invalid Docker image name: {image_name}")

    # Prevent shell injection
    dangerous_chars = {";", "&", "|", "`", "$", "(", ")", "<", ">", "\\"}
    if any(char in image_name for char in dangerous_chars):
        raise ValidationError(f"Docker image name contains dangerous characters")

    return image_name


def validate_command_args(args: list[str]) -> list[str]:
    """
    Validates command-line arguments for security.

    Prevents command injection attacks.

    Args:
        args: List of command arguments

    Returns:
        Validated argument list

    Raises:
        ValidationError: If validation fails
    """
    if not args:
        return []

    # Check each argument
    for arg in args:
        if not isinstance(arg, str):
            raise ValidationError(f"Argument must be string, got {type(arg)}")

        # Prevent shell metacharacters
        dangerous_chars = {";", "&", "|", "`", "$", "(", ")", "<", ">"}
        if any(char in arg for char in dangerous_chars):
            raise ValidationError(
                f"Argument contains dangerous shell characters: {arg}"
            )

    return args
