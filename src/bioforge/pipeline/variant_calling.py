"""
Variant Calling Module
Implements GATK best practices, FreeBayes, DeepVariant, and Strelka2
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import json
import shutil
import tempfile

logger = logging.getLogger(__name__)

@dataclass 
class VariantStats:
    """Statistics from variant calling"""
    total_variants: int = 0
    snps: int = 0
    indels: int = 0
    insertions: int = 0
    deletions: int = 0
    multi_allelic: int = 0
    homozygous: int = 0
    heterozygous: int = 0
    pass_filter: int = 0
    filtered: int = 0
    ti_tv_ratio: float = 0.0  # Transition/Transversion ratio
    mean_quality: float = 0.0
    mean_depth: float = 0.0


class VariantCaller:
    """
    Main class for variant calling operations
    """
    
    def __init__(self,
                 reference_fasta: str,
                 work_dir: str = "/tmp/variants",
                 threads: int = 4):
        """
        Initialize variant caller
        
        Args:
            reference_fasta: Path to reference FASTA
            work_dir: Working directory
            threads: Number of threads
        """
        self.reference_fasta = reference_fasta
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        
        # Check available variant callers
        self.callers = self._check_callers()
    
    def _check_callers(self) -> Dict[str, bool]:
        """Check which variant callers are available"""
        callers = {
            "gatk": shutil.which("gatk") is not None,
            "freebayes": shutil.which("freebayes") is not None,
            "bcftools": shutil.which("bcftools") is not None,
            "strelka": shutil.which("configureStrelkaGermlineWorkflow.py") is not None,
            "deepvariant": shutil.which("run_deepvariant") is not None,
            "vardict": shutil.which("vardict") is not None
        }
        
        for caller, available in callers.items():
            if available:
                logger.info(f"✓ {caller} available")
            else:
                logger.warning(f"✗ {caller} not found")
        
        return callers
    
    def _calculate_vcf_stats(self, vcf_file: str) -> VariantStats:
        """Calculate statistics from VCF file"""
        stats = VariantStats()
        
        try:
            transitions = 0
            transversions = 0
            total_quality = 0
            total_depth = 0
            
            with open(vcf_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) < 8:
                        continue
                    
                    stats.total_variants += 1
                    
                    # Parse variant info
                    ref = parts[3]
                    alt = parts[4]
                    qual = parts[5]
                    filter_field = parts[6]
                    info = parts[7]
                    
                    # Quality
                    if qual != '.':
                        total_quality += float(qual)
                    
                    # Filter status
                    if filter_field == 'PASS' or filter_field == '.':
                        stats.pass_filter += 1
                    else:
                        stats.filtered += 1
                    
                    # Variant type
                    if len(ref) == 1 and len(alt) == 1:
                        stats.snps += 1
                        # Ti/Tv calculation
                        if (ref in 'AG' and alt in 'AG') or (ref in 'CT' and alt in 'CT'):
                            transitions += 1
                        else:
                            transversions += 1
                    elif len(ref) > len(alt):
                        stats.deletions += 1
                        stats.indels += 1
                    elif len(ref) < len(alt):
                        stats.insertions += 1
                        stats.indels += 1
                    
                    # Parse INFO field for depth
                    for item in info.split(';'):
                        if item.startswith('DP='):
                            total_depth += int(item.split('=')[1])
                            break
                    
                    # Parse genotype if present
                    if len(parts) > 9:
                        gt = parts[9].split(':')[0]
                        if '/' in gt or '|' in gt:
                            alleles = gt.replace('|', '/').split('/')
                            if len(set(alleles)) == 1 and alleles[0] != '0':
                                stats.homozygous += 1
                            elif len(set(alleles)) > 1:
                                stats.heterozygous += 1
            
            # Calculate averages
            if stats.total_variants > 0:
                stats.mean_quality = total_quality / stats.total_variants
                stats.mean_depth = total_depth / stats.total_variants
            
            if transversions > 0:
                stats.ti_tv_ratio = transitions / transversions
                
        except Exception as e:
            logger.error(f"Failed to calculate VCF stats: {e}")
        
        return stats
    
    def _count_variants(self, vcf_file: str) -> int:
        """Count number of variants in VCF"""
        count = 0
        with open(vcf_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    count += 1
        return count


class JointGenotyping:
    """
    Joint genotyping for multiple samples
    """
    
    def __init__(self, reference_fasta: str, work_dir: str = "/tmp/joint"):
        self.reference_fasta = reference_fasta
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
    
    def create_genomicsdb(self,
                         gvcf_files: List[str],
                         interval: str,
                         db_path: str) -> bool:
        """
        Create GenomicsDB for joint calling
        
        Args:
            gvcf_files: List of GVCF files
            interval: Genomic interval
            db_path: Output database path
            
        Returns:
            Success status
        """
        try:
            # Create sample map
            sample_map = self.work_dir / "sample_map.txt"
            with open(sample_map, 'w') as f:
                for i, gvcf in enumerate(gvcf_files):
                    sample_name = f"sample_{i}"
                    f.write(f"{sample_name}\t{gvcf}\n")
            
            cmd = [
                "gatk", "GenomicsDBImport",
                "--genomicsdb-workspace-path", db_path,
                "--sample-name-map", str(sample_map),
                "--intervals", interval,
                "--reader-threads", "4"
            ]
            
            logger.info("Creating GenomicsDB")
            subprocess.run(cmd, check=True, capture_output=True)
            return True
            
        except Exception as e:
            logger.error(f"GenomicsDB creation failed: {e}")
            return False
