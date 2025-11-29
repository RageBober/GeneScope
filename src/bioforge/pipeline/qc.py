"""
Quality Control Module for Genomic Data
Handles FASTQ quality assessment, adapter trimming, and quality filtering
"""

import os
import subprocess
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass
import pandas as pd
import numpy as np
from datetime import datetime

logger = logging.getLogger(__name__)

@dataclass
class QCMetrics:
    """Quality control metrics for sequencing data"""
    total_reads: int = 0
    total_bases: int = 0
    q20_bases: int = 0  # Bases with quality >= 20
    q30_bases: int = 0  # Bases with quality >= 30
    gc_content: float = 0.0
    mean_quality: float = 0.0
    read_length_mean: float = 0.0
    read_length_std: float = 0.0
    adapter_content: float = 0.0
    duplication_rate: float = 0.0
    n_content: float = 0.0
    
    def to_dict(self) -> Dict:
        """Convert metrics to dictionary"""
        return {
            "total_reads": self.total_reads,
            "total_bases": self.total_bases,
            "q20_percentage": (self.q20_bases / self.total_bases * 100) if self.total_bases > 0 else 0,
            "q30_percentage": (self.q30_bases / self.total_bases * 100) if self.total_bases > 0 else 0,
            "gc_content": self.gc_content,
            "mean_quality": self.mean_quality,
            "read_length_mean": self.read_length_mean,
            "read_length_std": self.read_length_std,
            "adapter_content": self.adapter_content,
            "duplication_rate": self.duplication_rate,
            "n_content": self.n_content
        }
    
    def is_pass(self, thresholds: Dict = None) -> Tuple[bool, List[str]]:
        """
        Check if metrics pass quality thresholds
        
        Returns:
            Tuple of (pass/fail, list of failed criteria)
        """
        if thresholds is None:
            thresholds = {
                "min_reads": 1000000,  # 1M reads minimum
                "min_q30": 80,  # 80% Q30 bases
                "max_n_content": 5,  # Max 5% N bases
                "max_adapter": 10,  # Max 10% adapter
                "max_duplication": 50  # Max 50% duplication
            }
        
        failures = []
        
        if self.total_reads < thresholds["min_reads"]:
            failures.append(f"Low read count: {self.total_reads:,} < {thresholds['min_reads']:,}")
        
        q30_pct = (self.q30_bases / self.total_bases * 100) if self.total_bases > 0 else 0
        if q30_pct < thresholds["min_q30"]:
            failures.append(f"Low Q30: {q30_pct:.1f}% < {thresholds['min_q30']}%")
        
        if self.n_content > thresholds["max_n_content"]:
            failures.append(f"High N content: {self.n_content:.1f}% > {thresholds['max_n_content']}%")
        
        if self.adapter_content > thresholds["max_adapter"]:
            failures.append(f"High adapter content: {self.adapter_content:.1f}% > {thresholds['max_adapter']}%")
        
        if self.duplication_rate > thresholds["max_duplication"]:
            failures.append(f"High duplication: {self.duplication_rate:.1f}% > {thresholds['max_duplication']}%")
        
        return (len(failures) == 0, failures)


class QualityController:
    """
    Main class for quality control operations
    """
    
    def __init__(self, work_dir: str = "/tmp/qc", threads: int = 4):
        """
        Initialize QC controller
        
        Args:
            work_dir: Working directory for temporary files
            threads: Number of threads to use
        """
        self.work_dir = Path(work_dir)
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.threads = threads
        
        # Check for required tools
        self.tools = self._check_tools()
    
    def _check_tools(self) -> Dict[str, bool]:
        """Check availability of QC tools"""
        import shutil
        
        tools = {
            "fastqc": shutil.which("fastqc") is not None,
            "fastp": shutil.which("fastp") is not None,
            "trimmomatic": shutil.which("trimmomatic") is not None,
            "cutadapt": shutil.which("cutadapt") is not None,
            "multiqc": shutil.which("multiqc") is not None,
        }
        
        for tool, available in tools.items():
            if available:
                logger.info(f"✓ {tool} available")
            else:
                logger.warning(f"✗ {tool} not found")
        
        return tools
    
    def run_fastqc(self, 
                   fastq_files: List[str], 
                   output_dir: str) -> Dict[str, Any]:
        """
        Run FastQC quality assessment
        
        Args:
            fastq_files: List of FASTQ file paths
            output_dir: Output directory for FastQC reports
            
        Returns:
            Dict with QC results
        """
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        results = {
            "tool": "FastQC",
            "timestamp": datetime.now().isoformat(),
            "files": {}
        }
        
        if not self.tools.get("fastqc"):
            logger.warning("FastQC not available, using basic QC")
            return self._basic_qc(fastq_files)
        
        for fastq_file in fastq_files:
            try:
                cmd = [
                    "fastqc",
                    "-o", str(output_path),
                    "-t", str(self.threads),
                    "--quiet",
                    "--extract",
                    fastq_file
                ]
                
                logger.info(f"Running FastQC on {Path(fastq_file).name}")
                result = subprocess.run(cmd, capture_output=True, text=True)
                
                if result.returncode == 0:
                    # Parse FastQC results
                    base_name = Path(fastq_file).stem
                    if base_name.endswith('.fastq'):
                        base_name = base_name[:-6]
                    
                    data_file = output_path / f"{base_name}_fastqc" / "fastqc_data.txt"
                    if data_file.exists():
                        metrics = self._parse_fastqc_data(data_file)
                        results["files"][fastq_file] = metrics
                    else:
                        results["files"][fastq_file] = {"status": "completed", "data": "not parsed"}
                else:
                    logger.error(f"FastQC failed: {result.stderr}")
                    results["files"][fastq_file] = {"status": "failed", "error": result.stderr}
                    
            except Exception as e:
                logger.error(f"FastQC error for {fastq_file}: {e}")
                results["files"][fastq_file] = {"status": "error", "error": str(e)}
        
        return results
    
    def _parse_fastqc_data(self, data_file: Path) -> Dict:
        """Parse FastQC data file"""
        metrics = QCMetrics()
        current_module = None
        
        with open(data_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line.startswith(">>"):
                    parts = line[2:].split('\t')
                    current_module = parts[0]
                    continue
                
                if line.startswith("#") or not line:
                    continue
                
                # Parse basic statistics
                if current_module == "Basic Statistics":
                    if line.startswith("Total Sequences"):
                        metrics.total_reads = int(line.split('\t')[1])
                    elif line.startswith("Sequence length"):
                        length_str = line.split('\t')[1]
                        if '-' in length_str:
                            parts = length_str.split('-')
                            metrics.read_length_mean = (int(parts[0]) + int(parts[1])) / 2
                        else:
                            metrics.read_length_mean = float(length_str)
                    elif line.startswith("%GC"):
                        metrics.gc_content = float(line.split('\t')[1])
                
                # Parse per base sequence quality
                elif current_module == "Per base sequence quality":
                    parts = line.split('\t')
                    if len(parts) >= 2 and parts[0] != "Base":
                        try:
                            quality = float(parts[1])  # Mean quality
                            # This is simplified - real implementation would aggregate
                            metrics.mean_quality = quality
                        except:
                            pass
        
        return metrics.to_dict()
    
    def _basic_qc(self, fastq_files: List[str]) -> Dict:
        """Basic QC without external tools"""
        results = {
            "tool": "Basic QC",
            "timestamp": datetime.now().isoformat(),
            "files": {}
        }
        
        for fastq_file in fastq_files:
            metrics = QCMetrics()
            
            try:
                import gzip
                
                # Open file (handle gzipped)
                if fastq_file.endswith('.gz'):
                    f = gzip.open(fastq_file, 'rt')
                else:
                    f = open(fastq_file, 'r')
                
                line_count = 0
                total_quality = 0
                gc_count = 0
                n_count = 0
                
                for i, line in enumerate(f):
                    if i % 4 == 1:  # Sequence line
                        seq = line.strip()
                        metrics.total_bases += len(seq)
                        gc_count += seq.count('G') + seq.count('C')
                        n_count += seq.count('N')
                        metrics.total_reads += 1
                    elif i % 4 == 3:  # Quality line
                        qual = line.strip()
                        # Convert ASCII to Phred scores
                        for q in qual:
                            score = ord(q) - 33
                            total_quality += score
                            if score >= 20:
                                metrics.q20_bases += 1
                            if score >= 30:
                                metrics.q30_bases += 1
                
                f.close()
                
                # Calculate final metrics
                if metrics.total_bases > 0:
                    metrics.gc_content = (gc_count / metrics.total_bases) * 100
                    metrics.n_content = (n_count / metrics.total_bases) * 100
                    metrics.mean_quality = total_quality / metrics.total_bases
                
                results["files"][fastq_file] = metrics.to_dict()
                
            except Exception as e:
                logger.error(f"Basic QC failed for {fastq_file}: {e}")
                results["files"][fastq_file] = {"error": str(e)}
        
        return results
    
    def trim_adapters(self,
                     fastq_r1: str,
                     fastq_r2: Optional[str] = None,
                     output_dir: str = None,
                     adapter_r1: str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                     adapter_r2: str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT") -> Tuple[str, Optional[str], Dict]:
        """
        Trim adapters from FASTQ files using fastp or cutadapt
        
        Args:
            fastq_r1: Read 1 FASTQ file
            fastq_r2: Read 2 FASTQ file (optional)
            output_dir: Output directory
            adapter_r1: Adapter sequence for read 1
            adapter_r2: Adapter sequence for read 2
            
        Returns:
            Tuple of (trimmed_r1, trimmed_r2, metrics)
        """
        if output_dir is None:
            output_dir = str(self.work_dir / "trimmed")
        
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        base_name = Path(fastq_r1).stem.replace('.fastq', '').replace('.fq', '')
        
        if self.tools.get("fastp"):
            return self._trim_with_fastp(fastq_r1, fastq_r2, output_path, base_name)
        elif self.tools.get("cutadapt"):
            return self._trim_with_cutadapt(
                fastq_r1, fastq_r2, output_path, base_name, 
                adapter_r1, adapter_r2
            )
        else:
            logger.warning("No trimming tools available, skipping adapter trimming")
            return fastq_r1, fastq_r2, {"status": "skipped", "reason": "no tools"}
    
    def _trim_with_fastp(self, 
                        fastq_r1: str, 
                        fastq_r2: Optional[str],
                        output_path: Path,
                        base_name: str) -> Tuple[str, Optional[str], Dict]:
        """Trim using fastp - SECURE VERSION."""
        out_r1 = str(output_path / f"{base_name}_trimmed_R1.fastq.gz")
        out_r2 = str(output_path / f"{base_name}_trimmed_R2.fastq.gz") if fastq_r2 else None
        json_report = str(output_path / f"{base_name}_fastp.json")
        html_report = str(output_path / f"{base_name}_fastp.html")
        
        # Build secure command as list
        cmd = [
            "fastp",
            "-i", fastq_r1,
            "-o", out_r1,
            "--thread", str(self.threads),
            "--json", json_report,
            "--html", html_report,
            "--detect_adapter_for_pe",
            "--cut_front",
            "--cut_tail",
            "--cut_mean_quality", "20",
            "--length_required", "36",
            "--correction"
        ]
        
        if fastq_r2:
            cmd.extend(["-I", fastq_r2, "-O", out_r2])
        
        try:
            logger.info("Running fastp for adapter trimming")
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse fastp JSON report
            with open(json_report, 'r') as f:
                fastp_data = json.load(f)
            
            metrics = {
                "tool": "fastp",
                "status": "success",
                "before_filtering": {
                    "total_reads": fastp_data["summary"]["before_filtering"]["total_reads"],
                    "total_bases": fastp_data["summary"]["before_filtering"]["total_bases"],
                    "q30_rate": fastp_data["summary"]["before_filtering"]["q30_rate"]
                },
                "after_filtering": {
                    "total_reads": fastp_data["summary"]["after_filtering"]["total_reads"],
                    "total_bases": fastp_data["summary"]["after_filtering"]["total_bases"],
                    "q30_rate": fastp_data["summary"]["after_filtering"]["q30_rate"]
                },
                "filtering": {
                    "passed_filter_reads": fastp_data["filtering_result"]["passed_filter_reads"],
                    "low_quality_reads": fastp_data["filtering_result"]["low_quality_reads"],
                    "too_short_reads": fastp_data["filtering_result"]["too_many_N_reads"]
                }
            }
            
            return out_r1, out_r2, metrics
            
        except subprocess.CalledProcessError as e:
            logger.error(f"fastp failed with return code {e.returncode}")
            logger.error(f"fastp stderr: {e.stderr}")
            logger.error(f"fastp stdout: {e.stdout}")
            return fastq_r1, fastq_r2, {"status": "failed", "error": str(e)}
        except FileNotFoundError:
            logger.error("fastp JSON report not found after execution")
            return fastq_r1, fastq_r2, {"status": "failed", "error": "Report not generated"}
        except json.JSONDecodeError as e:
            logger.error(f"Failed to parse fastp JSON report: {e}")
            return fastq_r1, fastq_r2, {"status": "failed", "error": "Invalid JSON report"}
        except Exception as e:
            logger.error(f"Fastp trimming failed: {e}")
            return fastq_r1, fastq_r2, {"status": "failed", "error": str(e)}
    
    def _trim_with_cutadapt(self,
                           fastq_r1: str,
                           fastq_r2: Optional[str],
                           output_path: Path,
                           base_name: str,
                           adapter_r1: str,
                           adapter_r2: str) -> Tuple[str, Optional[str], Dict]:
        """Trim using cutadapt"""
        out_r1 = str(output_path / f"{base_name}_trimmed_R1.fastq.gz")
        
        cmd = [
            "cutadapt",
            "-a", adapter_r1,
            "-o", out_r1,
            "--minimum-length", "36",
            "--quality-cutoff", "20",
            "--cores", str(self.threads),
            fastq_r1
        ]
        
        if fastq_r2:
            out_r2 = str(output_path / f"{base_name}_trimmed_R2.fastq.gz")
            cmd = [
                "cutadapt",
                "-a", adapter_r1,
                "-A", adapter_r2,
                "-o", out_r1,
                "-p", out_r2,
                "--minimum-length", "36",
                "--quality-cutoff", "20",
                "--cores", str(self.threads),
                fastq_r1,
                fastq_r2
            ]
        else:
            out_r2 = None
        
        try:
            logger.info("Running cutadapt for adapter trimming")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse cutadapt output
                metrics = {
                    "tool": "cutadapt",
                    "status": "success",
                    "output": result.stdout
                }
                return out_r1, out_r2, metrics
            else:
                raise Exception(f"cutadapt failed: {result.stderr}")
                
        except Exception as e:
            logger.error(f"Cutadapt trimming failed: {e}")
            return fastq_r1, fastq_r2, {"status": "failed", "error": str(e)}
    
    def remove_duplicates(self, 
                         bam_file: str,
                         output_bam: str = None) -> Tuple[str, Dict]:
        """
        Remove PCR duplicates from BAM file
        
        Args:
            bam_file: Input BAM file
            output_bam: Output BAM file path
            
        Returns:
            Tuple of (deduplicated BAM path, metrics)
        """
        import shutil
        
        if output_bam is None:
            output_bam = bam_file.replace('.bam', '.dedup.bam')
        
        # Try different deduplication tools
        if shutil.which("picard"):
            return self._dedup_with_picard(bam_file, output_bam)
        elif shutil.which("samtools"):
            return self._dedup_with_samtools(bam_file, output_bam)
        else:
            logger.warning("No deduplication tools available")
            return bam_file, {"status": "skipped", "reason": "no tools"}
    
    def _dedup_with_picard(self, bam_file: str, output_bam: str) -> Tuple[str, Dict]:
        """Remove duplicates using Picard"""
        metrics_file = output_bam.replace('.bam', '.metrics.txt')
        
        cmd = [
            "picard", "MarkDuplicates",
            f"I={bam_file}",
            f"O={output_bam}",
            f"M={metrics_file}",
            "REMOVE_DUPLICATES=true",
            "CREATE_INDEX=true"
        ]
        
        try:
            logger.info("Running Picard MarkDuplicates")
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                # Parse metrics file
                metrics = {"tool": "picard", "status": "success"}
                if Path(metrics_file).exists():
                    with open(metrics_file, 'r') as f:
                        for line in f:
                            if line.startswith("LIBRARY"):
                                next(f)  # Skip header
                                data = next(f).strip().split('\t')
                                if len(data) >= 9:
                                    metrics["duplication_rate"] = float(data[8])
                                break
                
                return output_bam, metrics
            else:
                raise Exception(f"Picard failed: {result.stderr}")
                
        except Exception as e:
            logger.error(f"Picard deduplication failed: {e}")
            return bam_file, {"status": "failed", "error": str(e)}
    
    def _dedup_with_samtools(self, bam_file: str, output_bam: str) -> Tuple[str, Dict]:
        """Remove duplicates using samtools"""
        try:
            # Sort by name first
            name_sorted = output_bam.replace('.bam', '.namesort.bam')
            cmd1 = ["samtools", "sort", "-n", "-@", str(self.threads), "-o", name_sorted, bam_file]
            subprocess.run(cmd1, check=True)
            
            # Mark duplicates
            cmd2 = ["samtools", "fixmate", "-m", name_sorted, "-"]
            cmd3 = ["samtools", "sort", "-@", str(self.threads), "-"]
            cmd4 = ["samtools", "markdup", "-r", "-", output_bam]
            
            # Pipe commands
            p1 = subprocess.Popen(cmd2, stdout=subprocess.PIPE)
            p2 = subprocess.Popen(cmd3, stdin=p1.stdout, stdout=subprocess.PIPE)
            p3 = subprocess.Popen(cmd4, stdin=p2.stdout)
            p3.wait()
            
            # Index the output
            subprocess.run(["samtools", "index", output_bam], check=True)
            
            # Clean up
            Path(name_sorted).unlink(missing_ok=True)
            
            return output_bam, {"tool": "samtools", "status": "success"}
            
        except Exception as e:
            logger.error(f"Samtools deduplication failed: {e}")
            return bam_file, {"status": "failed", "error": str(e)}
    
    def generate_qc_report(self, 
                          qc_results: Dict,
                          output_file: str) -> str:
        """
        Generate comprehensive QC report
        
        Args:
            qc_results: QC results dictionary
            output_file: Output report file path
            
        Returns:
            Path to generated report
        """
        import json
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Generate HTML report
        html_content = self._generate_html_report(qc_results)
        
        if output_file.endswith('.json'):
            # Save as JSON
            with open(output_file, 'w') as f:
                json.dump(qc_results, f, indent=2, default=str)
        else:
            # Save as HTML
            with open(output_file, 'w') as f:
                f.write(html_content)
        
        logger.info(f"QC report generated: {output_file}")
        return output_file
    
    def _generate_html_report(self, qc_results: Dict) -> str:
        """Generate HTML QC report"""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>QC Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background: #4CAF50; color: white; padding: 20px; }}
                .metrics {{ margin: 20px 0; }}
                table {{ border-collapse: collapse; width: 100%; }}
                th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
                th {{ background: #f2f2f2; }}
                .pass {{ color: green; }}
                .fail {{ color: red; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>Quality Control Report</h1>
                <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="metrics">
                <h2>QC Summary</h2>
                <pre>{json.dumps(qc_results, indent=2, default=str)}</pre>
            </div>
        </body>
        </html>
        """
        return html
