"""
Pipeline Orchestrator for GenoScope

Manages the complete genomic analysis pipeline from raw reads to annotated variants.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Literal

from .qc import QualityController, QCMetrics
from .alignment import AlignmentEngine, AlignmentStats
from .variant_calling import VariantCaller, VariantStats
from .annotation import VariantAnnotator

logger = logging.getLogger(__name__)

PipelineStatus = Literal["pending", "running", "completed", "failed", "cancelled"]
AnalysisType = Literal["wgs", "wes", "panel", "rna-seq", "amplicon"]


@dataclass
class PipelineConfig:
    """Configuration for genomic analysis pipeline."""
    
    # Analysis type
    analysis_type: AnalysisType = "wgs"
    
    # Reference genome
    reference_genome: str = "hg38"
    reference_path: Optional[Path] = None
    
    # Tools to use
    alignment_tool: str = "bwa"
    variant_caller: str = "gatk"
    annotation_tool: str = "vep"
    
    # Quality thresholds
    min_base_quality: int = 20
    min_mapping_quality: int = 20
    min_variant_quality: float = 30.0
    min_depth: int = 10
    
    # Performance settings
    threads: int = 8
    memory_gb: int = 16
    
    # Output options
    keep_intermediate: bool = False
    generate_reports: bool = True
    
    # Additional parameters
    custom_params: Dict[str, Any] = field(default_factory=dict)


@dataclass
class PipelineResult:
    """Container for pipeline execution results."""
    
    pipeline_id: str
    status: PipelineStatus
    start_time: datetime
    end_time: Optional[datetime] = None
    
    # Input files
    input_files: List[Path] = field(default_factory=list)
    
    # Output files
    output_dir: Optional[Path] = None
    final_vcf: Optional[Path] = None
    final_bam: Optional[Path] = None
    qc_report: Optional[Path] = None
    
    # Metrics
    qc_metrics: Optional[QCMetrics] = None
    alignment_stats: Optional[AlignmentStats] = None
    variant_stats: Optional[VariantStats] = None
    
    # Errors and warnings
    errors: List[str] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    
    # Execution log
    log: List[Dict[str, Any]] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "pipeline_id": self.pipeline_id,
            "status": self.status,
            "start_time": self.start_time.isoformat(),
            "end_time": self.end_time.isoformat() if self.end_time else None,
            "duration_seconds": (self.end_time - self.start_time).total_seconds() if self.end_time else None,
            "input_files": [str(f) for f in self.input_files],
            "output_dir": str(self.output_dir) if self.output_dir else None,
            "final_vcf": str(self.final_vcf) if self.final_vcf else None,
            "final_bam": str(self.final_bam) if self.final_bam else None,
            "qc_report": str(self.qc_report) if self.qc_report else None,
            "qc_metrics": self.qc_metrics.to_dict() if self.qc_metrics else None,
            "alignment_stats": self.alignment_stats.to_dict() if self.alignment_stats else None,
            "variant_stats": self.variant_stats.to_dict() if self.variant_stats else None,
            "errors": self.errors,
            "warnings": self.warnings,
            "log": self.log
        }
    
    def save(self, output_file: Path) -> None:
        """Save results to JSON file."""
        with open(output_file, "w") as f:
            json.dump(self.to_dict(), f, indent=2)


class PipelineOrchestrator:
    """
    Orchestrates the complete genomic analysis pipeline.
    
    Coordinates QC, alignment, variant calling, and annotation steps.
    """
    
    def __init__(self,
                 config: PipelineConfig,
                 output_base_dir: Path = Path("./pipeline_outputs")):
        """
        Initialize pipeline orchestrator.
        
        Args:
            config: Pipeline configuration
            output_base_dir: Base directory for pipeline outputs
        """
        self.config = config
        self.output_base_dir = Path(output_base_dir)
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.qc_controller = QualityController(threads=config.threads)
        self.alignment_engine = AlignmentEngine(threads=config.threads)
        
        # Set reference path
        if config.reference_path:
            self.reference = config.reference_path
        else:
            # Use default reference location
            self.reference = self.alignment_engine.reference_dir / f"{config.reference_genome}.fa"
            if not self.reference.exists():
                raise FileNotFoundError(
                    f"Reference genome not found: {self.reference}. "
                    "Please download or specify path in config."
                )
        
        self.variant_caller = VariantCaller(
            reference=self.reference,
            threads=config.threads,
            java_mem=f"{config.memory_gb}G"
        )
        self.annotator = VariantAnnotator(threads=config.threads)
        
        # Current pipeline result
        self.result: Optional[PipelineResult] = None
    
    def run_pipeline(self,
                    fastq_r1: Path,
                    fastq_r2: Optional[Path] = None,
                    sample_name: str = "sample",
                    pipeline_id: Optional[str] = None) -> PipelineResult:
        """
        Run complete analysis pipeline.
        
        Args:
            fastq_r1: Read 1 FASTQ file
            fastq_r2: Read 2 FASTQ file (for paired-end)
            sample_name: Sample name
            pipeline_id: Unique pipeline ID
            
        Returns:
            Pipeline execution result
        """
        # Initialize result
        if pipeline_id is None:
            pipeline_id = f"{sample_name}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            
        self.result = PipelineResult(
            pipeline_id=pipeline_id,
            status="running",
            start_time=datetime.now(),
            input_files=[fastq_r1] + ([fastq_r2] if fastq_r2 else [])
        )
        
        # Create output directory
        output_dir = self.output_base_dir / pipeline_id
        output_dir.mkdir(parents=True, exist_ok=True)
        self.result.output_dir = output_dir
        
        try:
            # Step 1: Quality Control
            self._log("Starting quality control...")
            qc_result = self._run_qc(fastq_r1, fastq_r2, output_dir / "qc")
            
            # Step 2: Alignment
            self._log("Starting read alignment...")
            alignment_result = self._run_alignment(
                qc_result[0], qc_result[1], output_dir / "alignment", sample_name
            )
            
            # Step 3: Variant Calling
            self._log("Starting variant calling...")
            variant_result = self._run_variant_calling(
                alignment_result[0], output_dir / "variants"
            )
            
            # Step 4: Annotation
            self._log("Starting variant annotation...")
            annotation_result = self._run_annotation(
                variant_result[0], output_dir / "annotation"
            )
            
            # Step 5: Generate reports
            if self.config.generate_reports:
                self._log("Generating final reports...")
                self._generate_reports(output_dir / "reports")
            
            # Mark as completed
            self.result.status = "completed"
            self.result.end_time = datetime.now()
            self._log(f"Pipeline completed successfully in {(self.result.end_time - self.result.start_time).total_seconds():.1f} seconds")
            
        except Exception as e:
            logger.error(f"Pipeline failed: {str(e)}")
            self.result.status = "failed"
            self.result.errors.append(str(e))
            self.result.end_time = datetime.now()
            self._log(f"Pipeline failed: {str(e)}")
            
        finally:
            # Save results
            result_file = output_dir / "pipeline_result.json"
            self.result.save(result_file)
            
            # Clean up intermediate files if requested
            if not self.config.keep_intermediate:
                self._cleanup_intermediate_files(output_dir)
        
        return self.result
    
    def _run_qc(self,
               fastq_r1: Path,
               fastq_r2: Optional[Path],
               output_dir: Path) -> Tuple[Path, Optional[Path], QCMetrics]:
        """Run quality control step."""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run FastQC
        self.qc_controller.run_fastqc(
            [fastq_r1] + ([fastq_r2] if fastq_r2 else []),
            output_dir / "fastqc"
        )
        
        # Trim adapters and low-quality bases
        trimmed_r1, trimmed_r2, qc_metrics = self.qc_controller.trim_adapters(
            fastq_r1, fastq_r2, output_dir / "trimmed"
        )
        
        # Check QC metrics
        passes_qc, issues = qc_metrics.passes_qc()
        if not passes_qc:
            for issue in issues:
                self.result.warnings.append(f"QC warning: {issue}")
                self._log(f"QC warning: {issue}", level="warning")
        
        # Store metrics
        self.result.qc_metrics = qc_metrics
        
        # Generate QC report
        qc_report = output_dir / "qc_report.html"
        self.qc_controller.generate_qc_report(
            qc_metrics, qc_report, sample_name=self.result.pipeline_id
        )
        self.result.qc_report = qc_report
        
        self._log("Quality control completed")
        return trimmed_r1, trimmed_r2, qc_metrics
    
    def _run_alignment(self,
                      fastq_r1: Path,
                      fastq_r2: Optional[Path],
                      output_dir: Path,
                      sample_name: str) -> Tuple[Path, AlignmentStats]:
        """Run alignment step."""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Ensure reference is indexed
        self.alignment_engine.index_reference(self.reference, tool=self.config.alignment_tool)
        
        # Align reads
        bam_file, alignment_stats = self.alignment_engine.align_reads(
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            reference=self.reference,
            output_dir=output_dir,
            tool=self.config.alignment_tool,
            sample_name=sample_name
        )
        
        # Check alignment quality
        passes_qc, issues = alignment_stats.passes_qc()
        if not passes_qc:
            for issue in issues:
                self.result.warnings.append(f"Alignment warning: {issue}")
                self._log(f"Alignment warning: {issue}", level="warning")
        
        # Remove duplicates if requested
        if self.config.custom_params.get("remove_duplicates", True):
            bam_file, dup_rate = self.qc_controller.remove_duplicates(bam_file, output_dir)
            self._log(f"Removed duplicates: {dup_rate:.1f}% duplicate rate")
        
        # Store results
        self.result.final_bam = bam_file
        self.result.alignment_stats = alignment_stats
        
        self._log("Alignment completed")
        return bam_file, alignment_stats
    
    def _run_variant_calling(self,
                            bam_file: Path,
                            output_dir: Path) -> Tuple[Path, VariantStats]:
        """Run variant calling step."""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Call variants
        vcf_file, variant_stats = self.variant_caller.call_variants(
            bam_file=bam_file,
            output_dir=output_dir,
            caller=self.config.variant_caller,
            min_base_quality=self.config.min_base_quality,
            min_mapping_quality=self.config.min_mapping_quality
        )
        
        # Apply filters
        filtered_vcf = self.variant_caller.filter_variants(
            vcf_file=vcf_file,
            output_dir=output_dir,
            min_qual=self.config.min_variant_quality,
            min_depth=self.config.min_depth
        )
        
        # Check variant quality
        passes_qc, issues = variant_stats.quality_check()
        if not passes_qc:
            for issue in issues:
                self.result.warnings.append(f"Variant calling warning: {issue}")
                self._log(f"Variant calling warning: {issue}", level="warning")
        
        # Store results
        self.result.variant_stats = variant_stats
        
        self._log(f"Variant calling completed: {variant_stats.total_variants} variants found")
        return filtered_vcf, variant_stats
    
    def _run_annotation(self,
                       vcf_file: Path,
                       output_dir: Path) -> Path:
        """Run variant annotation step."""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Annotate variants
        annotated_vcf = self.annotator.annotate_vcf(
            vcf_file=vcf_file,
            output_dir=output_dir,
            tool=self.config.annotation_tool,
            genome_build=self.config.reference_genome
        )
        
        # Store final VCF
        self.result.final_vcf = annotated_vcf
        
        # Save annotation cache
        self.annotator.save_cache()
        
        self._log("Annotation completed")
        return annotated_vcf
    
    def _generate_reports(self, output_dir: Path) -> None:
        """Generate comprehensive analysis reports."""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate summary report
        summary_report = output_dir / "summary_report.html"
        self._generate_summary_report(summary_report)
        
        self._log("Reports generated")
    
    def _generate_summary_report(self, output_file: Path) -> None:
        """Generate HTML summary report."""
        html_template = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>GenoScope Analysis Report - {self.result.pipeline_id}</title>
            <style>
                body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #f5f5f5; }}
                .container {{ background: white; padding: 30px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
                h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
                h2 {{ color: #34495e; margin-top: 30px; }}
                .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin: 20px 0; }}
                .metric-card {{ 
                    background: #f8f9fa; 
                    padding: 20px; 
                    border-radius: 8px;
                    border-left: 4px solid #3498db;
                }}
                .metric-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}
                .metric-label {{ color: #7f8c8d; margin-top: 5px; }}
                .status-completed {{ color: #27ae60; font-weight: bold; }}
                .status-failed {{ color: #e74c3c; font-weight: bold; }}
                .warning {{ background: #fff3cd; border-left: 4px solid #ffc107; padding: 10px; margin: 10px 0; }}
                table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
                th {{ background: #34495e; color: white; padding: 12px; text-align: left; }}
                td {{ padding: 10px; border-bottom: 1px solid #ecf0f1; }}
                tr:hover {{ background: #f8f9fa; }}
            </style>
        </head>
        <body>
            <div class="container">
                <h1>🧬 GenoScope Analysis Report</h1>
                
                <h2>Pipeline Information</h2>
                <table>
                    <tr><td><strong>Pipeline ID:</strong></td><td>{self.result.pipeline_id}</td></tr>
                    <tr><td><strong>Status:</strong></td><td class="status-{self.result.status}">{self.result.status.upper()}</td></tr>
                    <tr><td><strong>Start Time:</strong></td><td>{self.result.start_time.strftime('%Y-%m-%d %H:%M:%S')}</td></tr>
                    <tr><td><strong>Duration:</strong></td><td>{(self.result.end_time - self.result.start_time).total_seconds():.1f} seconds</td></tr>
                    <tr><td><strong>Analysis Type:</strong></td><td>{self.config.analysis_type.upper()}</td></tr>
                    <tr><td><strong>Reference:</strong></td><td>{self.config.reference_genome}</td></tr>
                </table>
                
                <h2>Quality Control Metrics</h2>
                <div class="metric-grid">
                    <div class="metric-card">
                        <div class="metric-value">{self.result.qc_metrics.total_reads:,}</div>
                        <div class="metric-label">Total Reads</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.qc_metrics.mean_quality:.1f}</div>
                        <div class="metric-label">Mean Quality Score</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{(self.result.qc_metrics.q30_bases/self.result.qc_metrics.total_bases*100):.1f}%</div>
                        <div class="metric-label">Q30 Bases</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.qc_metrics.gc_content:.1f}%</div>
                        <div class="metric-label">GC Content</div>
                    </div>
                </div>
                
                <h2>Alignment Statistics</h2>
                <div class="metric-grid">
                    <div class="metric-card">
                        <div class="metric-value">{self.result.alignment_stats.mapping_rate:.1f}%</div>
                        <div class="metric-label">Mapping Rate</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.alignment_stats.mean_coverage:.1f}x</div>
                        <div class="metric-label">Mean Coverage</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.alignment_stats.properly_paired:,}</div>
                        <div class="metric-label">Properly Paired</div>
                    </div>
                </div>
                
                <h2>Variant Statistics</h2>
                <div class="metric-grid">
                    <div class="metric-card">
                        <div class="metric-value">{self.result.variant_stats.total_variants:,}</div>
                        <div class="metric-label">Total Variants</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.variant_stats.snps:,}</div>
                        <div class="metric-label">SNPs</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.variant_stats.indels:,}</div>
                        <div class="metric-label">Indels</div>
                    </div>
                    <div class="metric-card">
                        <div class="metric-value">{self.result.variant_stats.ti_tv_ratio:.2f}</div>
                        <div class="metric-label">Ti/Tv Ratio</div>
                    </div>
                </div>
                
                {"<h2>Warnings</h2>" + "".join([f'<div class="warning">{w}</div>' for w in self.result.warnings]) if self.result.warnings else ""}
                
                <h2>Output Files</h2>
                <table>
                    <tr><th>File Type</th><th>Path</th></tr>
                    <tr><td>Final VCF</td><td>{self.result.final_vcf}</td></tr>
                    <tr><td>Final BAM</td><td>{self.result.final_bam}</td></tr>
                    <tr><td>QC Report</td><td>{self.result.qc_report}</td></tr>
                </table>
                
                <p style="margin-top: 40px; color: #7f8c8d; text-align: center;">
                    <small>Report generated by GenoScope Pipeline v1.0</small>
                </p>
            </div>
        </body>
        </html>
        """
        
        output_file.write_text(html_template)
    
    def _cleanup_intermediate_files(self, output_dir: Path) -> None:
        """Clean up intermediate files to save space."""
        # Keep only final outputs
        # This is a placeholder - implement based on specific needs
        pass
    
    def _log(self, message: str, level: str = "info") -> None:
        """Add message to pipeline log."""
        log_entry = {
            "timestamp": datetime.now().isoformat(),
            "level": level,
            "message": message
        }
        
        if self.result:
            self.result.log.append(log_entry)
        
        if level == "error":
            logger.error(message)
        elif level == "warning":
            logger.warning(message)
        else:
            logger.info(message)
