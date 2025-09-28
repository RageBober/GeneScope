"""
Pipeline Validators Module for GenoScope

Provides validation functions for pipeline inputs and outputs.
"""

from __future__ import annotations

import hashlib
import logging
import os
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

logger = logging.getLogger(__name__)


class FileValidator:
    """Validates input files for genomic analysis."""
    
    # Supported file extensions
    VALID_EXTENSIONS = {
        ".fastq", ".fq", ".fastq.gz", ".fq.gz",  # FASTQ
        ".bam", ".sam", ".cram",  # Alignment
        ".vcf", ".vcf.gz", ".bcf",  # Variants
        ".bed", ".bedgraph",  # Regions
        ".gff", ".gff3", ".gtf",  # Annotations
        ".fa", ".fasta", ".fna"  # Reference
    }
    
    # Maximum file sizes (in bytes)
    MAX_SIZES = {
        ".fastq": 50 * 1024**3,  # 50 GB
        ".fastq.gz": 20 * 1024**3,  # 20 GB
        ".bam": 100 * 1024**3,  # 100 GB
        ".vcf": 10 * 1024**3,  # 10 GB
        ".vcf.gz": 5 * 1024**3,  # 5 GB
    }
    
    @staticmethod
    def validate_file_exists(file_path: Path) -> Tuple[bool, str]:
        """
        Check if file exists and is readable.
        
        Args:
            file_path: Path to file
            
        Returns:
            Tuple of (valid, message)
        """
        if not file_path.exists():
            return False, f"File does not exist: {file_path}"
            
        if not file_path.is_file():
            return False, f"Path is not a file: {file_path}"
            
        if not os.access(file_path, os.R_OK):
            return False, f"File is not readable: {file_path}"
            
        return True, "File is valid"
    
    @staticmethod
    def validate_file_format(file_path: Path) -> Tuple[bool, str]:
        """
        Validate file format based on extension.
        
        Args:
            file_path: Path to file
            
        Returns:
            Tuple of (valid, message)
        """
        # Get file extension(s)
        suffixes = "".join(file_path.suffixes)
        
        # Check if extension is valid
        valid_ext = False
        for ext in FileValidator.VALID_EXTENSIONS:
            if suffixes.endswith(ext):
                valid_ext = True
                break
                
        if not valid_ext:
            return False, f"Unsupported file format: {suffixes}"
            
        return True, f"Valid file format: {suffixes}"
    
    @staticmethod
    def validate_file_size(file_path: Path) -> Tuple[bool, str]:
        """
        Check if file size is within limits.
        
        Args:
            file_path: Path to file
            
        Returns:
            Tuple of (valid, message)
        """
        file_size = file_path.stat().st_size
        suffixes = "".join(file_path.suffixes)
        
        # Find matching max size
        max_size = None
        for ext, size in FileValidator.MAX_SIZES.items():
            if suffixes.endswith(ext):
                max_size = size
                break
                
        if max_size and file_size > max_size:
            return False, f"File too large: {file_size / 1024**3:.2f} GB (max: {max_size / 1024**3:.2f} GB)"
            
        return True, f"File size OK: {file_size / 1024**3:.2f} GB"
    
    @staticmethod
    def validate_fastq(file_path: Path, check_content: bool = False) -> Tuple[bool, str]:
        """
        Validate FASTQ file format.
        
        Args:
            file_path: Path to FASTQ file
            check_content: Whether to check file content
            
        Returns:
            Tuple of (valid, message)
        """
        if not check_content:
            return True, "FASTQ validation skipped"
            
        try:
            import gzip
            
            # Check if gzipped
            if file_path.suffix == ".gz":
                open_func = gzip.open
                mode = "rt"
            else:
                open_func = open
                mode = "r"
                
            # Read first few records
            with open_func(file_path, mode) as f:
                for i in range(4 * 10):  # Check first 10 records
                    line = f.readline()
                    if not line:
                        break
                        
                    if i % 4 == 0:  # Header line
                        if not line.startswith("@"):
                            return False, f"Invalid FASTQ header at line {i+1}"
                    elif i % 4 == 2:  # Quality header
                        if not line.startswith("+"):
                            return False, f"Invalid FASTQ quality header at line {i+1}"
                            
            return True, "Valid FASTQ file"
            
        except Exception as e:
            return False, f"Error reading FASTQ: {str(e)}"
    
    @staticmethod
    def calculate_md5(file_path: Path) -> str:
        """
        Calculate MD5 checksum of file.
        
        Args:
            file_path: Path to file
            
        Returns:
            MD5 checksum string
        """
        md5_hash = hashlib.md5()
        
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                md5_hash.update(chunk)
                
        return md5_hash.hexdigest()


class ParameterValidator:
    """Validates pipeline parameters."""
    
    @staticmethod
    def validate_threads(threads: int) -> Tuple[bool, str]:
        """
        Validate thread count.
        
        Args:
            threads: Number of threads
            
        Returns:
            Tuple of (valid, message)
        """
        if threads < 1:
            return False, "Thread count must be at least 1"
            
        max_threads = os.cpu_count() or 1
        if threads > max_threads * 2:  # Allow oversubscription up to 2x
            return False, f"Thread count ({threads}) exceeds available CPUs ({max_threads})"
            
        return True, f"Valid thread count: {threads}"
    
    @staticmethod
    def validate_memory(memory_gb: int) -> Tuple[bool, str]:
        """
        Validate memory allocation.
        
        Args:
            memory_gb: Memory in GB
            
        Returns:
            Tuple of (valid, message)
        """
        if memory_gb < 1:
            return False, "Memory must be at least 1 GB"
            
        # Get available memory (simplified)
        try:
            import psutil
            available_gb = psutil.virtual_memory().available / (1024**3)
            
            if memory_gb > available_gb:
                return False, f"Requested memory ({memory_gb} GB) exceeds available ({available_gb:.1f} GB)"
                
        except ImportError:
            logger.warning("psutil not installed, skipping memory check")
            
        return True, f"Valid memory allocation: {memory_gb} GB"
    
    @staticmethod
    def validate_quality_threshold(quality: int) -> Tuple[bool, str]:
        """
        Validate quality score threshold.
        
        Args:
            quality: Quality score
            
        Returns:
            Tuple of (valid, message)
        """
        if quality < 0:
            return False, "Quality score must be non-negative"
            
        if quality > 60:
            return False, "Quality score above 60 is unrealistic"
            
        if quality < 20:
            logger.warning(f"Low quality threshold: {quality}")
            
        return True, f"Valid quality threshold: {quality}"
    
    @staticmethod
    def validate_depth(depth: int) -> Tuple[bool, str]:
        """
        Validate read depth threshold.
        
        Args:
            depth: Read depth
            
        Returns:
            Tuple of (valid, message)
        """
        if depth < 1:
            return False, "Depth must be at least 1"
            
        if depth < 10:
            logger.warning(f"Low depth threshold: {depth}")
            
        return True, f"Valid depth threshold: {depth}"
    
    @staticmethod
    def validate_allele_frequency(freq: float) -> Tuple[bool, str]:
        """
        Validate allele frequency.
        
        Args:
            freq: Allele frequency
            
        Returns:
            Tuple of (valid, message)
        """
        if freq < 0 or freq > 1:
            return False, "Allele frequency must be between 0 and 1"
            
        return True, f"Valid allele frequency: {freq}"


class OutputValidator:
    """Validates pipeline outputs."""
    
    @staticmethod
    def validate_bam(bam_path: Path) -> Tuple[bool, str]:
        """
        Validate BAM file.
        
        Args:
            bam_path: Path to BAM file
            
        Returns:
            Tuple of (valid, message)
        """
        if not bam_path.exists():
            return False, f"BAM file does not exist: {bam_path}"
            
        # Check if indexed
        bai_path = bam_path.with_suffix(".bam.bai")
        if not bai_path.exists():
            bai_path = Path(str(bam_path) + ".bai")
            
        if not bai_path.exists():
            logger.warning(f"BAM index not found for {bam_path}")
            
        # Check with samtools if available
        try:
            import subprocess
            
            result = subprocess.run(
                ["samtools", "quickcheck", str(bam_path)],
                capture_output=True,
                text=True
            )
            
            if result.returncode != 0:
                return False, "BAM file is corrupted"
                
        except (FileNotFoundError, subprocess.SubprocessError):
            logger.debug("samtools not available for BAM validation")
            
        return True, "Valid BAM file"
    
    @staticmethod
    def validate_vcf(vcf_path: Path) -> Tuple[bool, str]:
        """
        Validate VCF file.
        
        Args:
            vcf_path: Path to VCF file
            
        Returns:
            Tuple of (valid, message)
        """
        if not vcf_path.exists():
            return False, f"VCF file does not exist: {vcf_path}"
            
        # Check header
        try:
            import gzip
            
            if vcf_path.suffix == ".gz":
                open_func = gzip.open
                mode = "rt"
            else:
                open_func = open
                mode = "r"
                
            with open_func(vcf_path, mode) as f:
                first_line = f.readline()
                if not first_line.startswith("##fileformat=VCF"):
                    return False, "Invalid VCF header"
                    
        except Exception as e:
            return False, f"Error reading VCF: {str(e)}"
            
        return True, "Valid VCF file"
    
    @staticmethod
    def validate_metrics(metrics: Dict[str, Any], thresholds: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """
        Validate pipeline metrics against thresholds.
        
        Args:
            metrics: Calculated metrics
            thresholds: Expected thresholds
            
        Returns:
            Tuple of (passes, failed_checks)
        """
        failed_checks = []
        
        # Check mapping rate
        if "mapping_rate" in metrics and "min_mapping_rate" in thresholds:
            if metrics["mapping_rate"] < thresholds["min_mapping_rate"]:
                failed_checks.append(
                    f"Low mapping rate: {metrics['mapping_rate']:.1f}% "
                    f"(threshold: {thresholds['min_mapping_rate']}%)"
                )
        
        # Check coverage
        if "mean_coverage" in metrics and "min_coverage" in thresholds:
            if metrics["mean_coverage"] < thresholds["min_coverage"]:
                failed_checks.append(
                    f"Low coverage: {metrics['mean_coverage']:.1f}x "
                    f"(threshold: {thresholds['min_coverage']}x)"
                )
        
        # Check variant count
        if "variant_count" in metrics:
            if "min_variants" in thresholds and metrics["variant_count"] < thresholds["min_variants"]:
                failed_checks.append(
                    f"Too few variants: {metrics['variant_count']} "
                    f"(expected: >{thresholds['min_variants']})"
                )
            
            if "max_variants" in thresholds and metrics["variant_count"] > thresholds["max_variants"]:
                failed_checks.append(
                    f"Too many variants: {metrics['variant_count']} "
                    f"(expected: <{thresholds['max_variants']})"
                )
        
        # Check Ti/Tv ratio
        if "ti_tv_ratio" in metrics:
            if metrics["ti_tv_ratio"] < 1.8 or metrics["ti_tv_ratio"] > 3.5:
                failed_checks.append(
                    f"Unusual Ti/Tv ratio: {metrics['ti_tv_ratio']:.2f} "
                    f"(expected: 1.8-3.5)"
                )
        
        return len(failed_checks) == 0, failed_checks


class PipelineValidator:
    """Main validator for complete pipeline."""
    
    def __init__(self):
        """Initialize pipeline validator."""
        self.file_validator = FileValidator()
        self.param_validator = ParameterValidator()
        self.output_validator = OutputValidator()
    
    def validate_input_files(self, 
                           fastq_r1: Path,
                           fastq_r2: Optional[Path] = None,
                           reference: Optional[Path] = None) -> Tuple[bool, List[str]]:
        """
        Validate all input files.
        
        Args:
            fastq_r1: Read 1 FASTQ
            fastq_r2: Read 2 FASTQ (optional)
            reference: Reference genome (optional)
            
        Returns:
            Tuple of (valid, error_messages)
        """
        errors = []
        
        # Validate Read 1
        valid, msg = self.file_validator.validate_file_exists(fastq_r1)
        if not valid:
            errors.append(msg)
        else:
            valid, msg = self.file_validator.validate_file_format(fastq_r1)
            if not valid:
                errors.append(msg)
                
            valid, msg = self.file_validator.validate_fastq(fastq_r1)
            if not valid:
                errors.append(msg)
        
        # Validate Read 2 if provided
        if fastq_r2:
            valid, msg = self.file_validator.validate_file_exists(fastq_r2)
            if not valid:
                errors.append(msg)
            else:
                valid, msg = self.file_validator.validate_file_format(fastq_r2)
                if not valid:
                    errors.append(msg)
                    
                valid, msg = self.file_validator.validate_fastq(fastq_r2)
                if not valid:
                    errors.append(msg)
        
        # Validate reference if provided
        if reference:
            valid, msg = self.file_validator.validate_file_exists(reference)
            if not valid:
                errors.append(msg)
            else:
                valid, msg = self.file_validator.validate_file_format(reference)
                if not valid:
                    errors.append(msg)
        
        return len(errors) == 0, errors
    
    def validate_parameters(self, config: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """
        Validate pipeline parameters.
        
        Args:
            config: Pipeline configuration
            
        Returns:
            Tuple of (valid, error_messages)
        """
        errors = []
        
        # Validate threads
        if "threads" in config:
            valid, msg = self.param_validator.validate_threads(config["threads"])
            if not valid:
                errors.append(msg)
        
        # Validate memory
        if "memory_gb" in config:
            valid, msg = self.param_validator.validate_memory(config["memory_gb"])
            if not valid:
                errors.append(msg)
        
        # Validate quality thresholds
        if "min_base_quality" in config:
            valid, msg = self.param_validator.validate_quality_threshold(config["min_base_quality"])
            if not valid:
                errors.append(msg)
        
        if "min_mapping_quality" in config:
            valid, msg = self.param_validator.validate_quality_threshold(config["min_mapping_quality"])
            if not valid:
                errors.append(msg)
        
        # Validate depth
        if "min_depth" in config:
            valid, msg = self.param_validator.validate_depth(config["min_depth"])
            if not valid:
                errors.append(msg)
        
        return len(errors) == 0, errors
    
    def validate_outputs(self,
                        bam_file: Optional[Path] = None,
                        vcf_file: Optional[Path] = None) -> Tuple[bool, List[str]]:
        """
        Validate pipeline outputs.
        
        Args:
            bam_file: Output BAM file
            vcf_file: Output VCF file
            
        Returns:
            Tuple of (valid, error_messages)
        """
        errors = []
        
        if bam_file:
            valid, msg = self.output_validator.validate_bam(bam_file)
            if not valid:
                errors.append(msg)
        
        if vcf_file:
            valid, msg = self.output_validator.validate_vcf(vcf_file)
            if not valid:
                errors.append(msg)
        
        return len(errors) == 0, errors
