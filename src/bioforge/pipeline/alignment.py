"""
Alignment Module for Genomic Data
Handles read alignment using BWA, Minimap2, STAR, and Bowtie2
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import json
import shutil

logger = logging.getLogger(__name__)

@dataclass
class AlignmentStats:
    """Statistics from alignment"""
    total_reads: int = 0
    mapped_reads: int = 0
    unmapped_reads: int = 0
    properly_paired: int = 0
    mapping_rate: float = 0.0
    average_mapq: float = 0.0
    duplicate_reads: int = 0
    secondary_alignments: int = 0
    supplementary_alignments: int = 0
    
    def calculate_rates(self):
        """Calculate mapping rates"""
        if self.total_reads > 0:
            self.mapping_rate = (self.mapped_reads / self.total_reads) * 100


class AlignmentEngine:
    """
    Main class for sequence alignment operations
    """
    
    # Reference genome configurations
    REFERENCE_GENOMES = {
        "GRCh38": {
            "fasta": "/data/references/GRCh38/GRCh38.fa",
            "bwa_index": "/data/references/GRCh38/GRCh38.fa",
            "minimap2_index": "/data/references/GRCh38/GRCh38.mmi",
            "star_index": "/data/references/GRCh38/star_index",
            "bowtie2_index": "/data/references/GRCh38/GRCh38",
            "gtf": "/data/references/GRCh38/gencode.v43.annotation.gtf"
        },
        "GRCh37": {
            "fasta": "/data/references/GRCh37/GRCh37.fa",
            "bwa_index": "/data/references/GRCh37/GRCh37.fa",
            "minimap2_index": "/data/references/GRCh37/GRCh37.mmi",
            "star_index": "/data/references/GRCh37/star_index",
            "bowtie2_index": "/data/references/GRCh37/GRCh37",
            "gtf": "/data/references/GRCh37/gencode.v19.annotation.gtf"
        }
    }
    
    def __init__(self, 
                 reference: str = "GRCh38",
                 work_dir: str = "/tmp/alignment",
                 threads: int = 4):
        """
        Initialize alignment engine
        
        Args:
            reference: Reference genome version
            work_dir: Working directory
            threads: Number of threads
        """
        self.reference = reference
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        
        # Get reference paths
        self.ref_config = self.REFERENCE_GENOMES.get(reference, {})
        
        # Check available aligners
        self.aligners = self._check_aligners()
    
    def _check_aligners(self) -> Dict[str, bool]:
        """Check which alignment tools are available"""
        aligners = {
            "bwa": shutil.which("bwa") is not None,
            "minimap2": shutil.which("minimap2") is not None,
            "star": shutil.which("STAR") is not None,
            "bowtie2": shutil.which("bowtie2") is not None,
            "samtools": shutil.which("samtools") is not None
        }
        
        for aligner, available in aligners.items():
            if available:
                logger.info(f"✓ {aligner} available")
            else:
                logger.warning(f"✗ {aligner} not found")
        
        return aligners
    
    def align_reads(self,
                   fastq_r1: str,
                   fastq_r2: Optional[str] = None,
                   output_bam: Optional[str] = None,
                   aligner: str = "bwa",
                   sample_name: str = "sample",
                   platform: str = "ILLUMINA",
                   read_type: str = "dna") -> Tuple[Optional[str], AlignmentStats]:
        """
        Align reads to reference genome
        
        Args:
            fastq_r1: Read 1 FASTQ file
            fastq_r2: Read 2 FASTQ file (optional)
            output_bam: Output BAM file path
            aligner: Aligner to use (bwa, minimap2, star, bowtie2)
            sample_name: Sample name for read group
            platform: Sequencing platform
            read_type: Type of reads (dna, rna, long)
            
        Returns:
            Tuple of (BAM file path, alignment statistics)
        """
        if output_bam is None:
            output_bam = str(self.work_dir / f"{sample_name}.bam")
        
        # Choose aligner based on read type
        if read_type == "rna" and self.aligners.get("star"):
            aligner = "star"
        elif read_type == "long" and self.aligners.get("minimap2"):
            aligner = "minimap2"
        
        # Check if aligner is available
        if not self.aligners.get(aligner):
            logger.error(f"{aligner} not available")
            # Try fallback aligners
            if self.aligners.get("bwa"):
                aligner = "bwa"
            elif self.aligners.get("minimap2"):
                aligner = "minimap2"
            else:
                return None, AlignmentStats()
        
        logger.info(f"Aligning with {aligner}")
        
        if aligner == "bwa":
            return self._align_bwa(fastq_r1, fastq_r2, output_bam, sample_name, platform)
        elif aligner == "minimap2":
            return self._align_minimap2(fastq_r1, fastq_r2, output_bam, sample_name, read_type)
        elif aligner == "star":
            return self._align_star(fastq_r1, fastq_r2, output_bam, sample_name)
        elif aligner == "bowtie2":
            return self._align_bowtie2(fastq_r1, fastq_r2, output_bam, sample_name)
        else:
            logger.error(f"Unknown aligner: {aligner}")
            return None, AlignmentStats()
    
    def _align_bwa(self,
                  fastq_r1: str,
                  fastq_r2: Optional[str],
                  output_bam: str,
                  sample_name: str,
                  platform: str) -> Tuple[Optional[str], AlignmentStats]:
        """Align using BWA-MEM"""
        ref_path = self.ref_config.get("bwa_index")
        if not ref_path or not Path(ref_path).exists():
            logger.error(f"BWA index not found: {ref_path}")
            return None, AlignmentStats()
        
        # Build read group
        read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:{platform}"
        
        # Build BWA command
        cmd = [
            "bwa", "mem",
            "-t", str(self.threads),
            "-R", read_group,
            ref_path,
            fastq_r1
        ]
        
        if fastq_r2:
            cmd.append(fastq_r2)
        
        try:
            # Run BWA and pipe to samtools
            sam_file = output_bam.replace('.bam', '.sam')
            
            logger.info(f"Running BWA alignment for {sample_name}")
            with open(sam_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                logger.error(f"BWA failed: {result.stderr}")
                return None, AlignmentStats()
            
            # Convert SAM to sorted BAM
            if self.aligners.get("samtools"):
                sorted_bam = self.sort_bam(sam_file, output_bam)
                Path(sam_file).unlink()  # Remove SAM file
                
                # Index BAM
                self.index_bam(sorted_bam)
                
                # Get alignment statistics
                stats = self.alignment_stats(sorted_bam)
                
                return sorted_bam, stats
            else:
                return sam_file, AlignmentStats()
                
        except Exception as e:
            logger.error(f"BWA alignment failed: {e}")
            return None, AlignmentStats()
    
    def _align_minimap2(self,
                       fastq_r1: str,
                       fastq_r2: Optional[str],
                       output_bam: str,
                       sample_name: str,
                       read_type: str) -> Tuple[Optional[str], AlignmentStats]:
        """Align using Minimap2"""
        ref_path = self.ref_config.get("minimap2_index")
        if not ref_path:
            ref_path = self.ref_config.get("fasta")
        
        if not ref_path or not Path(ref_path).exists():
            logger.error(f"Reference not found: {ref_path}")
            return None, AlignmentStats()
        
        # Choose preset based on read type
        preset_map = {
            "dna": "sr",  # Short reads
            "long": "map-ont",  # ONT long reads
            "pacbio": "map-pb",  # PacBio long reads
            "rna": "splice"  # RNA-seq
        }
        preset = preset_map.get(read_type, "sr")
        
        # Build minimap2 command
        cmd = [
            "minimap2",
            "-ax", preset,
            "-t", str(self.threads),
            "--MD",
            "-R", f"@RG\\tID:{sample_name}\\tSM:{sample_name}",
            ref_path,
            fastq_r1
        ]
        
        if fastq_r2:
            cmd.append(fastq_r2)
        
        try:
            sam_file = output_bam.replace('.bam', '.sam')
            
            logger.info(f"Running Minimap2 alignment for {sample_name}")
            with open(sam_file, 'w') as f:
                result = subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True)
            
            if result.returncode != 0:
                logger.error(f"Minimap2 failed: {result.stderr}")
                return None, AlignmentStats()
            
            # Convert and sort
            if self.aligners.get("samtools"):
                sorted_bam = self.sort_bam(sam_file, output_bam)
                Path(sam_file).unlink()
                self.index_bam(sorted_bam)
                stats = self.alignment_stats(sorted_bam)
                return sorted_bam, stats
            else:
                return sam_file, AlignmentStats()
                
        except Exception as e:
            logger.error(f"Minimap2 alignment failed: {e}")
            return None, AlignmentStats()
    
    def sort_bam(self, input_file: str, output_bam: str) -> str:
        """Sort BAM file using samtools"""
        cmd = [
            "samtools", "sort",
            "-@", str(self.threads),
            "-o", output_bam,
            input_file
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"BAM file sorted: {output_bam}")
            return output_bam
        except subprocess.CalledProcessError as e:
            logger.error(f"Sorting failed: {e.stderr}")
            return input_file
    
    def index_bam(self, bam_file: str) -> bool:
        """Index BAM file using samtools"""
        cmd = ["samtools", "index", bam_file]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"BAM file indexed: {bam_file}")
            return True
        except subprocess.CalledProcessError as e:
            logger.error(f"Indexing failed: {e.stderr}")
            return False
    
    def alignment_stats(self, bam_file: str) -> AlignmentStats:
        """Calculate alignment statistics using samtools flagstat"""
        stats = AlignmentStats()
        
        if not self.aligners.get("samtools"):
            return stats
        
        try:
            # Run samtools flagstat
            cmd = ["samtools", "flagstat", bam_file]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse flagstat output
            for line in result.stdout.split('\n'):
                if 'total' in line:
                    stats.total_reads = int(line.split()[0])
                elif 'mapped (' in line and 'secondary' not in line:
                    stats.mapped_reads = int(line.split()[0])
                elif 'properly paired' in line:
                    stats.properly_paired = int(line.split()[0])
                elif 'duplicates' in line:
                    stats.duplicate_reads = int(line.split()[0])
                elif 'secondary' in line:
                    stats.secondary_alignments = int(line.split()[0])
                elif 'supplementary' in line:
                    stats.supplementary_alignments = int(line.split()[0])
            
            stats.unmapped_reads = stats.total_reads - stats.mapped_reads
            stats.calculate_rates()
            
        except Exception as e:
            logger.error(f"Failed to get alignment stats: {e}")
        
        return stats
    
    def calculate_coverage(self, 
                          bam_file: str,
                          output_file: Optional[str] = None) -> Dict[str, Any]:
        """
        Calculate coverage statistics
        
        Args:
            bam_file: Input BAM file
            output_file: Output coverage file
            
        Returns:
            Coverage statistics
        """
        coverage_stats = {
            "mean_coverage": 0,
            "median_coverage": 0,
            "coverage_uniformity": 0,
            "percent_10x": 0,
            "percent_20x": 0,
            "percent_30x": 0
        }
        
        if output_file is None:
            output_file = bam_file.replace('.bam', '.coverage.txt')
        
        # Check for bedtools or mosdepth
        if shutil.which("mosdepth"):
            return self._coverage_mosdepth(bam_file, output_file)
        elif shutil.which("bedtools"):
            return self._coverage_bedtools(bam_file, output_file)
        elif self.aligners.get("samtools"):
            return self._coverage_samtools(bam_file, output_file)
        else:
            logger.warning("No coverage tools available")
            return coverage_stats
    
    def _coverage_samtools(self, bam_file: str, output_file: str) -> Dict[str, Any]:
        """Calculate coverage using samtools depth"""
        try:
            cmd = ["samtools", "depth", "-a", bam_file]
            
            with open(output_file, 'w') as f:
                subprocess.run(cmd, stdout=f, check=True)
            
            # Parse coverage file
            import numpy as np
            coverages = []
            
            with open(output_file, 'r') as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        coverages.append(int(parts[2]))
            
            if coverages:
                coverages = np.array(coverages)
                coverage_stats = {
                    "mean_coverage": float(np.mean(coverages)),
                    "median_coverage": float(np.median(coverages)),
                    "coverage_uniformity": float(np.std(coverages)),
                    "percent_10x": float(np.sum(coverages >= 10) / len(coverages) * 100),
                    "percent_20x": float(np.sum(coverages >= 20) / len(coverages) * 100),
                    "percent_30x": float(np.sum(coverages >= 30) / len(coverages) * 100)
                }
                return coverage_stats
            
        except Exception as e:
            logger.error(f"Coverage calculation failed: {e}")
        
        return {"error": "Coverage calculation failed"}
