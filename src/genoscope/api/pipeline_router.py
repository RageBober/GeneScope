"""
Pipeline API Router
Exposes genomics pipeline functionality via REST API
"""

from fastapi import APIRouter, HTTPException, BackgroundTasks, UploadFile, File, Form
from fastapi.responses import FileResponse
from typing import Optional, List, Dict, Any
from pathlib import Path
import uuid
import json
import logging
from datetime import datetime

from ..pipeline.main_pipeline import GenomicsPipeline, PipelineConfig, PipelineResult
from ..pipeline.qc import QualityController
from ..pipeline.alignment import AlignmentEngine
from ..pipeline.variant_calling import VariantCaller

logger = logging.getLogger(__name__)

# Create router
pipeline_router = APIRouter(prefix="/api/pipeline", tags=["pipeline"])

# Storage for pipeline jobs
PIPELINE_JOBS = {}
PIPELINE_DIR = Path("/data/pipeline_runs")
PIPELINE_DIR.mkdir(parents=True, exist_ok=True)

@pipeline_router.post("/submit")
async def submit_pipeline(
    background_tasks: BackgroundTasks,
    sample_name: str = Form(...),
    reference: str = Form("GRCh38"),
    aligner: str = Form("bwa"),
    variant_caller: str = Form("gatk"),
    fastq_r1: UploadFile = File(...),
    fastq_r2: Optional[UploadFile] = File(None),
    trim_adapters: bool = Form(True),
    remove_duplicates: bool = Form(True),
    annotate: bool = Form(True),
    threads: int = Form(4)
):
    """
    Submit a new pipeline job
    """
    # Generate job ID
    job_id = str(uuid.uuid4())
    job_dir = PIPELINE_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    
    try:
        # Save uploaded files
        input_dir = job_dir / "input"
        input_dir.mkdir(exist_ok=True)
        
        r1_path = input_dir / fastq_r1.filename
        with open(r1_path, 'wb') as f:
            content = await fastq_r1.read()
            f.write(content)
        
        r2_path = None
        if fastq_r2:
            r2_path = input_dir / fastq_r2.filename
            with open(r2_path, 'wb') as f:
                content = await fastq_r2.read()
                f.write(content)
        
        # Create pipeline configuration
        config = PipelineConfig(
            input_dir=str(input_dir),
            output_dir=str(job_dir / "output"),
            reference=reference,
            sample_name=sample_name,
            aligner=aligner,
            variant_caller=variant_caller,
            trim_adapters=trim_adapters,
            remove_duplicates=remove_duplicates,
            annotate=annotate,
            threads=threads
        )
        
        # Store job info
        PIPELINE_JOBS[job_id] = {
            "status": "submitted",
            "submitted_at": datetime.now().isoformat(),
            "config": config.__dict__,
            "result": None
        }
        
        # Run pipeline in background
        background_tasks.add_task(
            run_pipeline_task,
            job_id,
            config,
            str(r1_path),
            str(r2_path) if r2_path else None
        )
        
        return {
            "job_id": job_id,
            "status": "submitted",
            "message": f"Pipeline job submitted for {sample_name}"
        }
        
    except Exception as e:
        logger.error(f"Failed to submit pipeline: {e}")
        raise HTTPException(status_code=500, detail=str(e))

def run_pipeline_task(job_id: str, config: PipelineConfig, r1_path: str, r2_path: Optional[str]):
    """Run pipeline in background"""
    try:
        PIPELINE_JOBS[job_id]["status"] = "running"
        PIPELINE_JOBS[job_id]["started_at"] = datetime.now().isoformat()
        
        # Initialize and run pipeline
        pipeline = GenomicsPipeline(config)
        result = pipeline.run(r1_path, r2_path)
        
        # Store result
        PIPELINE_JOBS[job_id]["status"] = result.status
        PIPELINE_JOBS[job_id]["result"] = result.__dict__
        PIPELINE_JOBS[job_id]["completed_at"] = datetime.now().isoformat()
        
        # Save result to file
        result_file = Path(config.output_dir) / "pipeline_result.json"
        with open(result_file, 'w') as f:
            json.dump(result.__dict__, f, indent=2, default=str)
            
    except Exception as e:
        logger.error(f"Pipeline failed for job {job_id}: {e}")
        PIPELINE_JOBS[job_id]["status"] = "failed"
        PIPELINE_JOBS[job_id]["error"] = str(e)

@pipeline_router.get("/status/{job_id}")
async def get_pipeline_status(job_id: str):
    """Get pipeline job status"""
    if job_id not in PIPELINE_JOBS:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = PIPELINE_JOBS[job_id]
    
    return {
        "job_id": job_id,
        "status": job["status"],
        "submitted_at": job.get("submitted_at"),
        "started_at": job.get("started_at"),
        "completed_at": job.get("completed_at"),
        "sample_name": job["config"]["sample_name"],
        "error": job.get("error")
    }

@pipeline_router.get("/result/{job_id}")
async def get_pipeline_result(job_id: str):
    """Get pipeline job results"""
    if job_id not in PIPELINE_JOBS:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = PIPELINE_JOBS[job_id]
    
    if job["status"] not in ["completed", "failed_qc", "failed_alignment", "failed_variant_calling", "failed"]:
        return {
            "job_id": job_id,
            "status": job["status"],
            "message": "Pipeline still running"
        }
    
    return {
        "job_id": job_id,
        "status": job["status"],
        "result": job.get("result"),
        "error": job.get("error")
    }

@pipeline_router.get("/download/{job_id}/{file_type}")
async def download_result_file(job_id: str, file_type: str):
    """
    Download result files
    
    file_type options:
    - vcf: Variant call file
    - bam: Aligned BAM file
    - qc_report: QC report
    - full_report: Complete pipeline report
    """
    if job_id not in PIPELINE_JOBS:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = PIPELINE_JOBS[job_id]
    result = job.get("result")
    
    if not result:
        raise HTTPException(status_code=400, detail="Results not available yet")
    
    job_dir = PIPELINE_DIR / job_id / "output"
    
    file_map = {
        "vcf": result.get("vcf_file"),
        "annotated_vcf": result.get("annotated_vcf"),
        "bam": result.get("aligned_bam"),
        "qc_report": str(job_dir / "reports" / "qc_report.html"),
        "full_report": str(job_dir / "reports" / "pipeline_result.json")
    }
    
    file_path = file_map.get(file_type)
    
    if not file_path or not Path(file_path).exists():
        raise HTTPException(status_code=404, detail=f"File type {file_type} not found")
    
    return FileResponse(
        file_path,
        filename=Path(file_path).name,
        media_type="application/octet-stream"
    )

@pipeline_router.get("/jobs")
async def list_pipeline_jobs(
    status: Optional[str] = None,
    limit: int = 10
):
    """List pipeline jobs"""
    jobs = []
    
    for job_id, job in PIPELINE_JOBS.items():
        if status and job["status"] != status:
            continue
        
        jobs.append({
            "job_id": job_id,
            "status": job["status"],
            "sample_name": job["config"]["sample_name"],
            "submitted_at": job.get("submitted_at"),
            "completed_at": job.get("completed_at")
        })
    
    # Sort by submission time
    jobs.sort(key=lambda x: x.get("submitted_at", ""), reverse=True)
    
    return {
        "total": len(jobs),
        "jobs": jobs[:limit]
    }

@pipeline_router.post("/qc-only")
async def run_qc_only(
    fastq_file: UploadFile = File(...),
    background_tasks: BackgroundTasks = None
):
    """Run QC analysis only"""
    try:
        # Save uploaded file
        temp_dir = Path("/tmp") / str(uuid.uuid4())
        temp_dir.mkdir(parents=True, exist_ok=True)
        
        file_path = temp_dir / fastq_file.filename
        with open(file_path, 'wb') as f:
            content = await fastq_file.read()
            f.write(content)
        
        # Run QC
        qc = QualityController(work_dir=str(temp_dir))
        results = qc.run_fastqc([str(file_path)], str(temp_dir))
        
        # Parse metrics
        metrics = {}
        for file, data in results.get("files", {}).items():
            if isinstance(data, dict):
                metrics = data
                break
        
        return {
            "filename": fastq_file.filename,
            "metrics": metrics,
            "status": "completed"
        }
        
    except Exception as e:
        logger.error(f"QC analysis failed: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@pipeline_router.get("/reference-genomes")
async def list_reference_genomes():
    """List available reference genomes"""
    from ..pipeline.alignment import AlignmentEngine
    
    references = []
    for ref_name, ref_config in AlignmentEngine.REFERENCE_GENOMES.items():
        # Check if reference files exist
        fasta_exists = Path(ref_config.get("fasta", "")).exists()
        
        references.append({
            "name": ref_name,
            "available": fasta_exists,
            "paths": ref_config
        })
    
    return {"references": references}

@pipeline_router.get("/tools")
async def check_available_tools():
    """Check which bioinformatics tools are available"""
    import shutil
    
    tools = {
        # QC tools
        "fastqc": shutil.which("fastqc") is not None,
        "fastp": shutil.which("fastp") is not None,
        "multiqc": shutil.which("multiqc") is not None,
        
        # Aligners
        "bwa": shutil.which("bwa") is not None,
        "minimap2": shutil.which("minimap2") is not None,
        "star": shutil.which("STAR") is not None,
        "bowtie2": shutil.which("bowtie2") is not None,
        
        # Variant callers
        "gatk": shutil.which("gatk") is not None,
        "freebayes": shutil.which("freebayes") is not None,
        "bcftools": shutil.which("bcftools") is not None,
        "strelka": shutil.which("configureStrelkaGermlineWorkflow.py") is not None,
        
        # Utilities
        "samtools": shutil.which("samtools") is not None,
        "bedtools": shutil.which("bedtools") is not None,
        "picard": shutil.which("picard") is not None,
        
        # Annotation
        "vep": shutil.which("vep") is not None,
        "snpeff": shutil.which("snpEff") is not None,
        "annovar": shutil.which("annotate_variation.pl") is not None,
    }
    
    available = [tool for tool, present in tools.items() if present]
    missing = [tool for tool, present in tools.items() if not present]
    
    return {
        "available": available,
        "missing": missing,
        "total_available": len(available),
        "total_missing": len(missing),
        "details": tools
    }
