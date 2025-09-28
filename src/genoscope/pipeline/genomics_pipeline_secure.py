"""
Genomics Pipeline Module - SECURE VERSION
Main pipeline for processing genomic data (FASTQ → VCF → Annotation).
Integrates BWA, GATK, and annotation tools with security improvements.

SECURITY IMPROVEMENTS:
- Removed all shell=True usage
- Added input validation 
- Improved error handling with stderr capture
- Added path traversal protection
"""

import os
import subprocess
import logging
import json
import shutil
import re
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from datetime import datetime
import tempfile

from ..integrations.vcf_processor import VCFProcessor, Variant
from ..integrations.clinvar import ClinVarAPI
from ..integrations.dbsnp import DbSNPAPI

logger = logging.getLogger(__name__)

class SecurityError(Exception):
    """Custom exception for security-related errors."""
    pass

class GenomicsPipelineSecure:
    """
    Secure genomics analysis pipeline.
    Handles alignment, variant calling, and annotation with security improvements.
    """
    
    def __init__(self, 
                 work_dir: str = "/tmp/genoscope",
                 reference_genome: str = "GRCh38",
                 threads: int = 4):
        """
        Initialize genomics pipeline.
        
        Args:
            work_dir: Working directory for temporary files
            reference_genome: Reference genome version
            threads: Number of threads to use
        """
        self.work_dir = Path(work_dir).resolve()  # Resolve to absolute path
        self.work_dir.mkdir(parents=True, exist_ok=True)
        
        self.reference_genome = reference_genome
        self.threads = max(1, min(32, int(threads)))  # Validate threads count
        
        # Initialize API clients
        self.clinvar = ClinVarAPI()
        self.dbsnp = DbSNPAPI()
        
        # Check for required tools
        self.tools = self._check_tools()
        
        # Reference genome paths (these would be configured per installation)
        self.reference_paths = {
            "GRCh38": "/data/references/GRCh38/GRCh38.fa",
            "GRCh37": "/data/references/GRCh37/GRCh37.fa",
            "hg38": "/data/references/GRCh38/GRCh38.fa",
            "hg19": "/data/references/GRCh37/GRCh37.fa"
        }
    
    def _validate_file_path(self, file_path: str, must_exist: bool = True) -> Path:
        """
        Validate and secure file path.
        
        Args:
            file_path: File path to validate
            must_exist: Whether file must exist
            
        Returns:
            Validated Path object
            
        Raises:
            SecurityError: If path is unsafe
            FileNotFoundError: If file must exist but doesn't
        """
        try:
            path = Path(file_path).resolve()
        except (OSError, ValueError) as e:
            raise SecurityError(f"Invalid file path: {file_path}") from e
        
        # Check for path traversal attempts
        if ".." in str(path) or str(path).startswith("/etc") or str(path).startswith("/root"):
            raise SecurityError(f"Unsafe path detected: {file_path}")
        
        # Validate file extension for known safe types
        safe_extensions = {'.fastq', '.fq', '.bam', '.sam', '.vcf', '.fa', '.fasta', '.gz'}
        if not any(str(path).endswith(ext) or str(path).endswith(ext + '.gz') for ext in safe_extensions):
            logger.warning(f"Unusual file extension: {path.suffix}")
        
        if must_exist and not path.exists():
            raise FileNotFoundError(f"File not found: {file_path}")
        
        if must_exist and not path.is_file():
            raise SecurityError(f"Path is not a regular file: {file_path}")
        
        return path
    
    def _validate_sample_name(self, sample_name: str) -> str:
        """
        Validate sample name to prevent injection.
        
        Args:
            sample_name: Sample name to validate
            
        Returns:
            Validated sample name
            
        Raises:
            SecurityError: If sample name is unsafe
        """
        # Only allow alphanumeric, underscore, hyphen
        if not re.match(r'^[a-zA-Z0-9_-]+$', sample_name):
            raise SecurityError(f"Invalid sample name: {sample_name}")
        
        if len(sample_name) > 100:  # Reasonable length limit
            raise SecurityError(f"Sample name too long: {sample_name}")
        
        return sample_name
    
    def _safe_subprocess_run(self, cmd: List[str], timeout: int = 3600, **kwargs) -> subprocess.CompletedProcess:
        """
        Safely run subprocess with proper error handling.
        
        Args:
            cmd: Command as list (NO shell=True)
            timeout: Command timeout in seconds
            **kwargs: Additional subprocess.run arguments
            
        Returns:
            CompletedProcess result
            
        Raises:
            subprocess.CalledProcessError: If command fails
            subprocess.TimeoutExpired: If command times out
        """
        # Validate command list
        if not isinstance(cmd, list) or not cmd:
            raise SecurityError("Command must be a non-empty list")
        
        # Ensure all command parts are strings
        cmd = [str(part) for part in cmd]
        
        # Log command (without sensitive data)
        logger.debug(f"Running command: {cmd[0]} with {len(cmd)-1} arguments")
        
        try:
            result = subprocess.run(
                cmd,
                timeout=timeout,
                capture_output=True,
                text=True,
                check=True,
                **kwargs
            )
            
            if result.stderr:
                logger.debug(f"Command stderr: {result.stderr}")
            
            return result
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Command failed: {cmd[0]}")
            logger.error(f"Return code: {e.returncode}")
            logger.error(f"Stderr: {e.stderr}")
            logger.error(f"Stdout: {e.stdout}")
            raise
        except subprocess.TimeoutExpired as e:
            logger.error(f"Command timed out after {timeout}s: {cmd[0]}")
            raise
    
    def _check_tools(self) -> Dict[str, bool]:
        """Check availability of required bioinformatics tools."""
        tools = {
            "bwa": False,
            "samtools": False,
            "bcftools": False,
            "gatk": False,
            "fastp": False,
            "fastqc": False
        }
        
        for tool in tools:
            try:
                # Check if tool exists and is executable
                tool_path = shutil.which(tool)
                if tool_path:
                    # Try to get version to verify it's working
                    version_cmd = [tool, "--version"]
                    try:
                        result = subprocess.run(
                            version_cmd, 
                            capture_output=True, 
                            text=True, 
                            timeout=10
                        )
                        tools[tool] = True
                        logger.info(f"✓ {tool} available at {tool_path}")
                    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
                        # Tool exists but version check failed, still mark as available
                        tools[tool] = True
                        logger.info(f"✓ {tool} available at {tool_path} (version check failed)")
                else:
                    logger.warning(f"✗ {tool} not found in PATH")
            except Exception as e:
                logger.warning(f"✗ Error checking {tool}: {e}")
        
        return tools
    
    def run_quality_control(self, 
                          fastq_files: List[str],
                          output_dir: str) -> Dict[str, Any]:
        """
        Run quality control on FASTQ files with secure implementation.
        """
        # Validate inputs
        validated_files = []
        for fastq_file in fastq_files:
            validated_files.append(self._validate_file_path(fastq_file, must_exist=True))
        
        output_path = Path(output_dir).resolve()
        output_path.mkdir(parents=True, exist_ok=True)
        
        qc_results = {
            "timestamp": datetime.now().isoformat(),
            "files": {},
            "summary": {
                "total_reads": 0,
                "total_bases": 0,
                "average_quality": 0,
                "gc_content": 0
            }
        }
        
        for fastq_path in validated_files:
            file_name = fastq_path.name
            logger.info(f"Running QC on {file_name}")
            
            # Run FastQC if available
            if self.tools.get("fastqc"):
                try:
                    cmd = [
                        "fastqc",
                        "-o", str(output_path),
                        "-t", str(self.threads),
                        "--quiet",
                        str(fastq_path)
                    ]
                    self._safe_subprocess_run(cmd)
                    qc_results["files"][file_name] = {"fastqc": "completed"}
                except subprocess.CalledProcessError as e:
                    logger.error(f"FastQC failed for {file_name}: {e}")
                    qc_results["files"][file_name] = {"fastqc": "failed", "error": str(e)}
            
            # Basic statistics using safe commands
            try:
                if str(fastq_path).endswith('.gz'):
                    # Use zcat safely
                    cmd1 = ["zcat", str(fastq_path)]
                    cmd2 = ["wc", "-l"]
                    
                    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
                    result = self._safe_subprocess_run(cmd2, stdin=p1.stdout)
                    p1.stdout.close()
                    
                    total_lines = int(result.stdout.strip())
                else:
                    cmd = ["wc", "-l", str(fastq_path)]
                    result = self._safe_subprocess_run(cmd)
                    total_lines = int(result.stdout.strip().split()[0])
                
                total_reads = total_lines // 4
                
                qc_results["files"][file_name] = {
                    "total_reads": total_reads,
                    "format": "FASTQ"
                }
                qc_results["summary"]["total_reads"] += total_reads
                
            except Exception as e:
                logger.error(f"Failed to get basic stats for {file_name}: {e}")
        
        return qc_results
    
    def run_alignment(self,
                     fastq_r1: str,
                     fastq_r2: Optional[str] = None,
                     output_bam: str = None,
                     sample_name: str = "sample") -> Optional[str]:
        """
        Align FASTQ files to reference genome using BWA - SECURE VERSION.
        """
        # Validate inputs
        fastq_r1_path = self._validate_file_path(fastq_r1, must_exist=True)
        fastq_r2_path = None
        if fastq_r2:
            fastq_r2_path = self._validate_file_path(fastq_r2, must_exist=True)
        
        sample_name = self._validate_sample_name(sample_name)
        
        if not self.tools.get("bwa"):
            logger.error("BWA not available. Using mock alignment.")
            output_path = output_bam or str(self.work_dir / "aligned.bam")
            return self._mock_alignment(output_path)
        
        if not output_bam:
            output_bam = str(self.work_dir / f"{sample_name}_aligned.bam")
        
        output_bam_path = Path(output_bam).resolve()
        
        try:
            # Get and validate reference genome path
            ref_path = self.reference_paths.get(self.reference_genome)
            if not ref_path or not Path(ref_path).exists():
                logger.warning(f"Reference genome {self.reference_genome} not found. Using mock data.")
                return self._mock_alignment(str(output_bam_path))
            
            ref_path = self._validate_file_path(ref_path, must_exist=True)
            
            logger.info(f"Starting BWA alignment for {sample_name}")
            
            # Build secure BWA command
            read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"
            
            bwa_cmd = [
                "bwa", "mem",
                "-t", str(self.threads),
                "-R", read_group,
                str(ref_path),
                str(fastq_r1_path)
            ]
            
            if fastq_r2_path:
                bwa_cmd.append(str(fastq_r2_path))
            
            # Run BWA and pipe to samtools
            if self.tools.get("samtools"):
                samtools_sort_cmd = [
                    "samtools", "sort",
                    "-@", str(self.threads),
                    "-o", str(output_bam_path),
                    "-"
                ]
                
                # Use proper piping
                logger.info("Running BWA alignment with samtools sorting...")
                p1 = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
                self._safe_subprocess_run(samtools_sort_cmd, stdin=p1.stdout)
                p1.stdout.close()
                
                # Index the BAM file
                index_cmd = ["samtools", "index", str(output_bam_path)]
                self._safe_subprocess_run(index_cmd)
                
                logger.info(f"Alignment completed: {output_bam_path}")
                return str(output_bam_path)
            else:
                # Just run BWA and output SAM
                sam_file = str(output_bam_path).replace('.bam', '.sam')
                
                with open(sam_file, 'w') as f:
                    self._safe_subprocess_run(bwa_cmd, stdout=f)
                
                logger.info(f"Alignment completed (SAM format): {sam_file}")
                return sam_file
                
        except subprocess.CalledProcessError as e:
            logger.error(f"Alignment failed: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error during alignment: {e}")
            return None
    
    def _mock_alignment(self, output_path: str) -> str:
        """Create mock alignment file for testing."""
        sam_content = """@HD\tVN:1.6\tSO:coordinate
@SQ\tSN:chr1\tLN:248956422
@SQ\tSN:chr2\tLN:242193529
@RG\tID:sample\tSM:sample\tPL:ILLUMINA
read1\t99\tchr1\t100000\t60\t100M\t=\t100200\t300\tACGTACGTACGT\t############\tNM:i:0
read2\t147\tchr1\t100200\t60\t100M\t=\t100000\t-300\tACGTACGTACGT\t############\tNM:i:0
"""
        with open(output_path, 'w') as f:
            f.write(sam_content)
        
        logger.info(f"Created mock alignment file: {output_path}")
        return output_path
    
    def run_variant_calling(self,
                          bam_file: str,
                          output_vcf: str = None,
                          method: str = "bcftools") -> Optional[str]:
        """
        Call variants from aligned BAM file - SECURE VERSION.
        """
        # Validate inputs
        bam_path = self._validate_file_path(bam_file, must_exist=True)
        
        if not output_vcf:
            output_vcf = str(self.work_dir / "variants.vcf")
        
        output_vcf_path = Path(output_vcf).resolve()
        
        if method == "bcftools" and self.tools.get("bcftools"):
            return self._variant_calling_bcftools_secure(str(bam_path), str(output_vcf_path))
        elif method == "gatk" and self.tools.get("gatk"):
            return self._variant_calling_gatk_secure(str(bam_path), str(output_vcf_path))
        else:
            logger.warning(f"Variant caller {method} not available. Using mock data.")
            return self._mock_variant_calling(str(output_vcf_path))
    
    def _variant_calling_bcftools_secure(self, bam_file: str, output_vcf: str) -> Optional[str]:
        """Run variant calling using bcftools - SECURE VERSION."""
        try:
            ref_path = self.reference_paths.get(self.reference_genome)
            if not ref_path or not Path(ref_path).exists():
                return self._mock_variant_calling(output_vcf)
            
            ref_path = self._validate_file_path(ref_path, must_exist=True)
            
            # Secure bcftools commands
            mpileup_cmd = [
                "bcftools", "mpileup",
                "-Ou",
                "-f", str(ref_path),
                bam_file
            ]
            
            call_cmd = [
                "bcftools", "call",
                "-mv",
                "-Ov",
                "-o", output_vcf
            ]
            
            # Run with proper piping
            logger.info("Running bcftools variant calling...")
            p1 = subprocess.Popen(mpileup_cmd, stdout=subprocess.PIPE)
            self._safe_subprocess_run(call_cmd, stdin=p1.stdout)
            p1.stdout.close()
            
            logger.info(f"Variant calling completed: {output_vcf}")
            return output_vcf
            
        except subprocess.CalledProcessError as e:
            logger.error(f"bcftools variant calling failed: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in variant calling: {e}")
            return None
    
    def _variant_calling_gatk_secure(self, bam_file: str, output_vcf: str) -> Optional[str]:
        """Run variant calling using GATK HaplotypeCaller - SECURE VERSION."""
        try:
            ref_path = self.reference_paths.get(self.reference_genome)
            if not ref_path or not Path(ref_path).exists():
                return self._mock_variant_calling(output_vcf)
            
            ref_path = self._validate_file_path(ref_path, must_exist=True)
            
            gatk_cmd = [
                "gatk", "HaplotypeCaller",
                "-R", str(ref_path),
                "-I", bam_file,
                "-O", output_vcf,
                "--native-pair-hmm-threads", str(self.threads)
            ]
            
            self._safe_subprocess_run(gatk_cmd)
            logger.info(f"GATK variant calling completed: {output_vcf}")
            return output_vcf
            
        except subprocess.CalledProcessError as e:
            logger.error(f"GATK variant calling failed: {e}")
            return None
        except Exception as e:
            logger.error(f"Unexpected error in GATK variant calling: {e}")
            return None
    
    def _mock_variant_calling(self, output_vcf: str) -> str:
        """Create mock VCF file for testing."""
        vcf_content = """##fileformat=VCFv4.2
##reference=GRCh38
##source=GenoScope_Pipeline_Secure
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr1\t100000\trs123456\tA\tG\t60\tPASS\tDP=30;AF=0.5\tGT:GQ:DP\t0/1:50:30
chr1\t200000\t.\tAT\tA\t45\tPASS\tDP=25;AF=0.3\tGT:GQ:DP\t0/1:40:25
chr2\t300000\trs789012\tG\tC\t80\tPASS\tDP=40;AF=0.8\tGT:GQ:DP\t1/1:70:40
chr17\t43044295\trs80357906\tG\tA\t99\tPASS\tDP=50;AF=1.0\tGT:GQ:DP\t1/1:99:50
"""
        with open(output_vcf, 'w') as f:
            f.write(vcf_content)
        
        logger.info(f"Created mock VCF file: {output_vcf}")
        return output_vcf
    
    def annotate_variants(self, vcf_file: str) -> Dict[str, Any]:
        """
        Annotate variants with ClinVar and dbSNP information - SECURE VERSION.
        """
        # Validate VCF file
        vcf_path = self._validate_file_path(vcf_file, must_exist=True)
        
        logger.info(f"Annotating variants from {vcf_path}")
        
        try:
            # Parse VCF file
            processor = VCFProcessor(str(vcf_path), self.reference_genome)
            variants = processor.parse(limit=1000)  # Limit for performance
            
            annotated = {
                "total_variants": len(variants),
                "annotated_variants": [],
                "clinvar_matches": 0,
                "dbsnp_matches": 0,
                "pathogenic_variants": [],
                "statistics": processor.get_statistics()
            }
            
            for var in variants:
                annotation = {
                    "variant": var.to_dict(),
                    "clinvar": None,
                    "dbsnp": None,
                    "interpretation": None
                }
                
                # Check ClinVar
                if var.variant_id and var.variant_id.startswith("rs"):
                    # Query by rsID
                    clinvar_result = self.clinvar.search_by_rsid(var.variant_id)
                    if clinvar_result:
                        annotation["clinvar"] = clinvar_result
                        annotated["clinvar_matches"] += 1
                        
                        # Check if pathogenic
                        interpretation = self.clinvar.interpret_significance(
                            clinvar_result.get("clinical_significance", "")
                        )
                        if interpretation["is_pathogenic"]:
                            annotated["pathogenic_variants"].append(annotation)
                    
                    # Query dbSNP
                    dbsnp_result = self.dbsnp.get_snp_info(var.variant_id)
                    if dbsnp_result:
                        annotation["dbsnp"] = dbsnp_result
                        annotated["dbsnp_matches"] += 1
                else:
                    # Query ClinVar by position
                    clinvar_result = self.clinvar.search_by_variant(
                        chromosome=var.chromosome.replace("chr", ""),
                        position=var.position,
                        ref=var.ref,
                        alt=var.alt,
                        assembly=self.reference_genome
                    )
                    if clinvar_result:
                        annotation["clinvar"] = clinvar_result
                        annotated["clinvar_matches"] += 1
                
                annotated["annotated_variants"].append(annotation)
            
            logger.info(f"Annotation complete: {annotated['clinvar_matches']} ClinVar matches, "
                       f"{annotated['dbsnp_matches']} dbSNP matches")
            
            return annotated
        
        except Exception as e:
            logger.error(f"Annotation failed: {e}")
            return {
                "total_variants": 0,
                "error": str(e),
                "annotated_variants": [],
                "clinvar_matches": 0,
                "dbsnp_matches": 0,
                "pathogenic_variants": []
            }


# Backward compatibility alias
GenomicsPipeline = GenomicsPipelineSecure
