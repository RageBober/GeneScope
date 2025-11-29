"""ARQ Worker Tasks for BioForge.

Real genomics pipeline tasks that run in the background worker process.
Uses subprocess to call external bioinformatics tools (BWA, GATK, samtools).
"""

import asyncio
import logging
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bioforge.config import settings
from bioforge.core.exceptions import PipelineError
from bioforge.db import get_session_context
from bioforge.models import Job, JobStatus, Sample

logger = logging.getLogger("bioforge.workers")


# ─────────────────────────────────────────────────────────────────────────────
# Utility Functions
# ─────────────────────────────────────────────────────────────────────────────


def check_tool(tool: str) -> bool:
    """Check if a bioinformatics tool is available."""
    return shutil.which(tool) is not None


async def run_command(
    cmd: list[str],
    cwd: Path | None = None,
    timeout: int = 3600,
) -> tuple[int, str, str]:
    """Run a shell command asynchronously."""
    logger.info(f"Running: {' '.join(cmd)}")

    process = await asyncio.create_subprocess_exec(
        *cmd,
        cwd=cwd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE,
    )

    try:
        stdout, stderr = await asyncio.wait_for(
            process.communicate(),
            timeout=timeout,
        )
        return (
            process.returncode or 0,
            stdout.decode() if stdout else "",
            stderr.decode() if stderr else "",
        )
    except asyncio.TimeoutError:
        process.kill()
        raise PipelineError(f"Command timed out after {timeout}s: {' '.join(cmd)}")


async def update_job_progress(
    session: Any,
    job: Job,
    step: str,
    progress: int,
) -> None:
    """Update job progress in database."""
    job.current_step = step
    job.progress = progress
    await session.commit()
    logger.info(f"Job {job.id}: {step} ({progress}%)")


# ─────────────────────────────────────────────────────────────────────────────
# Pipeline Tasks
# ─────────────────────────────────────────────────────────────────────────────


async def run_alignment(ctx: dict, job_id: int) -> dict[str, Any]:
    """
    Run BWA-MEM2 alignment pipeline for a sample.

    Steps:
    1. Index reference (if needed)
    2. Align reads with BWA-MEM2
    3. Sort and convert to BAM
    4. Mark duplicates
    5. Index final BAM

    Args:
        ctx: ARQ context (contains redis connection)
        job_id: ID of the job to execute

    Returns:
        Result dict with status and output path
    """
    logger.info(f"Starting alignment job {job_id}")

    async with get_session_context() as session:
        from sqlalchemy import select

        # Get job and sample
        result = await session.execute(
            select(Job, Sample).join(Sample).where(Job.id == job_id)
        )
        row = result.first()
        if not row:
            logger.error(f"Job {job_id} not found")
            return {"status": "error", "message": "Job not found"}

        job, sample = row[0], row[1]

        # Update status to running
        job.status = JobStatus.RUNNING.value
        job.started_at = datetime.now(timezone.utc)
        await update_job_progress(session, job, "Initializing", 0)

        try:
            # Setup paths
            output_dir = (
                settings.storage.results_dir
                / str(sample.project_id)
                / str(sample.id)
                / str(job.id)
            )
            output_dir.mkdir(parents=True, exist_ok=True)

            # Get input files from sample
            fastq_r1 = Path(sample.fastq_r1) if sample.fastq_r1 else None
            fastq_r2 = Path(sample.fastq_r2) if sample.fastq_r2 else None

            if not fastq_r1 or not fastq_r1.exists():
                raise PipelineError(f"FASTQ R1 not found: {fastq_r1}")

            # Reference genome
            reference = settings.storage.references_dir / "GRCh38" / "GRCh38.fa"
            if not reference.exists():
                # Use test reference for development
                reference = settings.storage.references_dir / "test_reference.fasta"

            # Check tools
            use_bwa_mem2 = check_tool("bwa-mem2")
            use_bwa = check_tool("bwa")
            use_samtools = check_tool("samtools")

            if not (use_bwa_mem2 or use_bwa):
                raise PipelineError("Neither bwa-mem2 nor bwa found in PATH")
            if not use_samtools:
                raise PipelineError("samtools not found in PATH")

            bwa_cmd = "bwa-mem2" if use_bwa_mem2 else "bwa"
            threads = str(settings.pipeline.threads)
            sample_name = sample.name or f"sample_{sample.id}"

            # ─────────────────────────────────────────────────────────────────
            # Step 1: Index reference (if needed)
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Checking reference index", 5)

            index_suffix = ".bwt" if not use_bwa_mem2 else ".bwt.2bit.64"
            if not Path(str(reference) + index_suffix).exists():
                await update_job_progress(session, job, "Indexing reference", 10)
                returncode, stdout, stderr = await run_command(
                    [bwa_cmd, "index", str(reference)],
                    timeout=7200,
                )
                if returncode != 0:
                    raise PipelineError(f"BWA index failed: {stderr}")

            # ─────────────────────────────────────────────────────────────────
            # Step 2: Align reads
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Aligning reads", 20)

            sam_file = output_dir / f"{sample_name}.sam"
            read_group = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

            align_cmd = [
                bwa_cmd,
                "mem",
                "-t",
                threads,
                "-R",
                read_group,
                str(reference),
                str(fastq_r1),
            ]
            if fastq_r2 and fastq_r2.exists():
                align_cmd.append(str(fastq_r2))

            returncode, stdout, stderr = await run_command(align_cmd, timeout=7200)
            if returncode != 0:
                raise PipelineError(f"BWA alignment failed: {stderr}")

            # Write SAM output
            sam_file.write_text(stdout)
            await update_job_progress(session, job, "Alignment complete", 50)

            # ─────────────────────────────────────────────────────────────────
            # Step 3: Sort and convert to BAM
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Sorting BAM", 60)

            sorted_bam = output_dir / f"{sample_name}.sorted.bam"
            returncode, stdout, stderr = await run_command(
                [
                    "samtools",
                    "sort",
                    "-@",
                    threads,
                    "-o",
                    str(sorted_bam),
                    str(sam_file),
                ],
                timeout=3600,
            )
            if returncode != 0:
                raise PipelineError(f"Samtools sort failed: {stderr}")

            # Remove SAM file
            sam_file.unlink(missing_ok=True)

            # ─────────────────────────────────────────────────────────────────
            # Step 4: Mark duplicates (using samtools)
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Marking duplicates", 75)

            dedup_bam = output_dir / f"{sample_name}.sorted.dedup.bam"

            # Use samtools markdup (GATK MarkDuplicates alternative)
            if check_tool("samtools"):
                # First fixmate for markdup
                fixmate_bam = output_dir / f"{sample_name}.fixmate.bam"
                returncode, _, stderr = await run_command(
                    [
                        "samtools",
                        "fixmate",
                        "-m",
                        str(sorted_bam),
                        str(fixmate_bam),
                    ],
                    timeout=1800,
                )

                if returncode == 0:
                    # Sort by coordinate again
                    sorted_fixmate = output_dir / f"{sample_name}.sorted.fixmate.bam"
                    await run_command(
                        [
                            "samtools",
                            "sort",
                            "-o",
                            str(sorted_fixmate),
                            str(fixmate_bam),
                        ]
                    )
                    fixmate_bam.unlink(missing_ok=True)

                    # Mark duplicates
                    returncode, _, stderr = await run_command(
                        [
                            "samtools",
                            "markdup",
                            str(sorted_fixmate),
                            str(dedup_bam),
                        ],
                        timeout=1800,
                    )
                    sorted_fixmate.unlink(missing_ok=True)

                if returncode != 0:
                    # Fallback: just use sorted BAM
                    logger.warning(f"markdup failed, using sorted BAM: {stderr}")
                    shutil.copy(sorted_bam, dedup_bam)
            else:
                shutil.copy(sorted_bam, dedup_bam)

            sorted_bam.unlink(missing_ok=True)

            # ─────────────────────────────────────────────────────────────────
            # Step 5: Index BAM
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Indexing BAM", 90)

            returncode, _, stderr = await run_command(
                ["samtools", "index", str(dedup_bam)],
                timeout=600,
            )
            if returncode != 0:
                logger.warning(f"BAM indexing failed: {stderr}")

            # ─────────────────────────────────────────────────────────────────
            # Step 6: Calculate stats
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Calculating statistics", 95)

            stats = {}
            returncode, stdout, _ = await run_command(
                ["samtools", "flagstat", str(dedup_bam)]
            )
            if returncode == 0:
                stats["flagstat"] = stdout

            # Update sample with BAM path
            sample.bam_path = str(dedup_bam)

            # Mark job complete
            job.status = JobStatus.COMPLETED.value
            job.progress = 100
            job.current_step = "Completed"
            job.completed_at = datetime.now(timezone.utc)
            job.result_path = str(output_dir)
            job.result_data = {"bam": str(dedup_bam), "stats": stats}
            await session.commit()

            logger.info(f"Alignment job {job_id} completed: {dedup_bam}")
            return {"status": "completed", "output": str(dedup_bam), "stats": stats}

        except Exception as e:
            logger.exception(f"Alignment job {job_id} failed: {e}")
            job.status = JobStatus.FAILED.value
            job.error_message = str(e)
            job.completed_at = datetime.now(timezone.utc)
            await session.commit()
            return {"status": "failed", "error": str(e)}


async def run_variant_calling(ctx: dict, job_id: int) -> dict[str, Any]:
    """
    Run variant calling pipeline using bcftools or GATK.

    Steps:
    1. Generate pileup
    2. Call variants
    3. Filter variants
    4. Index VCF

    Args:
        ctx: ARQ context
        job_id: ID of the job to execute

    Returns:
        Result dict with status and output path
    """
    logger.info(f"Starting variant calling job {job_id}")

    async with get_session_context() as session:
        from sqlalchemy import select

        result = await session.execute(
            select(Job, Sample).join(Sample).where(Job.id == job_id)
        )
        row = result.first()
        if not row:
            return {"status": "error", "message": "Job not found"}

        job, sample = row[0], row[1]

        job.status = JobStatus.RUNNING.value
        job.started_at = datetime.now(timezone.utc)
        await update_job_progress(session, job, "Initializing", 0)

        try:
            # Setup paths
            output_dir = (
                settings.storage.results_dir
                / str(sample.project_id)
                / str(sample.id)
                / str(job.id)
            )
            output_dir.mkdir(parents=True, exist_ok=True)

            # Get BAM file
            bam_path = Path(sample.bam_path) if sample.bam_path else None
            if not bam_path or not bam_path.exists():
                raise PipelineError(f"BAM file not found: {bam_path}")

            # Reference genome
            reference = settings.storage.references_dir / "GRCh38" / "GRCh38.fa"
            if not reference.exists():
                reference = settings.storage.references_dir / "test_reference.fasta"

            sample_name = sample.name or f"sample_{sample.id}"

            # Check available tools
            use_gatk = check_tool("gatk")
            use_bcftools = check_tool("bcftools")

            if not (use_gatk or use_bcftools):
                raise PipelineError("Neither GATK nor bcftools found in PATH")

            # ─────────────────────────────────────────────────────────────────
            # Variant Calling with bcftools (faster) or GATK
            # ─────────────────────────────────────────────────────────────────

            raw_vcf = output_dir / f"{sample_name}.raw.vcf"
            filtered_vcf = output_dir / f"{sample_name}.filtered.vcf.gz"

            if use_bcftools:
                # bcftools mpileup + call pipeline
                await update_job_progress(session, job, "Generating pileup", 20)

                pileup_cmd = [
                    "bcftools",
                    "mpileup",
                    "-f",
                    str(reference),
                    "-Ou",
                    str(bam_path),
                ]

                call_cmd = [
                    "bcftools",
                    "call",
                    "-mv",
                    "-Ov",
                    "-o",
                    str(raw_vcf),
                ]

                # Run mpileup | call
                await update_job_progress(session, job, "Calling variants", 40)

                mpileup_proc = await asyncio.create_subprocess_exec(
                    *pileup_cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )

                call_proc = await asyncio.create_subprocess_exec(
                    *call_cmd,
                    stdin=mpileup_proc.stdout,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                )

                await mpileup_proc.wait()
                await call_proc.wait()

                if call_proc.returncode != 0:
                    _, stderr = await call_proc.communicate()
                    raise PipelineError(f"bcftools call failed: {stderr.decode()}")

                # Filter variants
                await update_job_progress(session, job, "Filtering variants", 70)

                returncode, _, stderr = await run_command(
                    [
                        "bcftools",
                        "filter",
                        "-s",
                        "LowQual",
                        "-e",
                        "QUAL<20 || DP<5",
                        "-Oz",
                        "-o",
                        str(filtered_vcf),
                        str(raw_vcf),
                    ]
                )
                if returncode != 0:
                    logger.warning(f"bcftools filter failed: {stderr}")
                    # Compress raw VCF as fallback
                    await run_command(
                        ["bgzip", "-c", str(raw_vcf)],
                    )

            elif use_gatk:
                # GATK HaplotypeCaller
                await update_job_progress(session, job, "Running HaplotypeCaller", 30)

                returncode, _, stderr = await run_command(
                    [
                        "gatk",
                        "HaplotypeCaller",
                        "-R",
                        str(reference),
                        "-I",
                        str(bam_path),
                        "-O",
                        str(raw_vcf),
                        "--native-pair-hmm-threads",
                        str(settings.pipeline.threads),
                    ],
                    timeout=7200,
                )
                if returncode != 0:
                    raise PipelineError(f"HaplotypeCaller failed: {stderr}")

                # Filter with GATK
                await update_job_progress(session, job, "Filtering variants", 70)

                returncode, _, stderr = await run_command(
                    [
                        "gatk",
                        "VariantFiltration",
                        "-R",
                        str(reference),
                        "-V",
                        str(raw_vcf),
                        "-O",
                        str(filtered_vcf),
                        "--filter-expression",
                        "QUAL < 30.0",
                        "--filter-name",
                        "LowQual",
                        "--filter-expression",
                        "DP < 10",
                        "--filter-name",
                        "LowDepth",
                    ]
                )

            # ─────────────────────────────────────────────────────────────────
            # Index VCF
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Indexing VCF", 90)

            if check_tool("tabix"):
                await run_command(["tabix", "-p", "vcf", str(filtered_vcf)])
            elif check_tool("bcftools"):
                await run_command(["bcftools", "index", str(filtered_vcf)])

            # Clean up
            raw_vcf.unlink(missing_ok=True)

            # ─────────────────────────────────────────────────────────────────
            # Calculate stats
            # ─────────────────────────────────────────────────────────────────
            await update_job_progress(session, job, "Calculating statistics", 95)

            stats = {}
            if check_tool("bcftools"):
                returncode, stdout, _ = await run_command(
                    ["bcftools", "stats", str(filtered_vcf)]
                )
                if returncode == 0:
                    # Parse key stats
                    for line in stdout.split("\n"):
                        if line.startswith("SN"):
                            parts = line.split("\t")
                            if len(parts) >= 4:
                                key = parts[2].rstrip(":")
                                value = parts[3]
                                stats[key] = value

            # Update sample
            sample.vcf_path = str(filtered_vcf)

            # Mark job complete
            job.status = JobStatus.COMPLETED.value
            job.progress = 100
            job.current_step = "Completed"
            job.completed_at = datetime.now(timezone.utc)
            job.result_path = str(output_dir)
            job.result_data = {"vcf": str(filtered_vcf), "stats": stats}
            await session.commit()

            logger.info(f"Variant calling job {job_id} completed: {filtered_vcf}")
            return {
                "status": "completed",
                "output": str(filtered_vcf),
                "stats": stats,
            }

        except Exception as e:
            logger.exception(f"Variant calling job {job_id} failed")
            job.status = JobStatus.FAILED.value
            job.error_message = str(e)
            job.completed_at = datetime.now(timezone.utc)
            await session.commit()
            return {"status": "failed", "error": str(e)}


async def run_qc(ctx: dict, job_id: int) -> dict[str, Any]:
    """
    Run quality control on FASTQ files.

    Uses fastp for QC and trimming.

    Args:
        ctx: ARQ context
        job_id: ID of the job to execute

    Returns:
        Result dict with QC metrics
    """
    logger.info(f"Starting QC job {job_id}")

    async with get_session_context() as session:
        from sqlalchemy import select

        result = await session.execute(
            select(Job, Sample).join(Sample).where(Job.id == job_id)
        )
        row = result.first()
        if not row:
            return {"status": "error", "message": "Job not found"}

        job, sample = row[0], row[1]

        job.status = JobStatus.RUNNING.value
        job.started_at = datetime.now(timezone.utc)
        await update_job_progress(session, job, "Initializing", 0)

        try:
            output_dir = (
                settings.storage.results_dir
                / str(sample.project_id)
                / str(sample.id)
                / str(job.id)
            )
            output_dir.mkdir(parents=True, exist_ok=True)

            fastq_r1 = Path(sample.fastq_r1) if sample.fastq_r1 else None
            fastq_r2 = Path(sample.fastq_r2) if sample.fastq_r2 else None

            if not fastq_r1 or not fastq_r1.exists():
                raise PipelineError(f"FASTQ R1 not found: {fastq_r1}")

            sample_name = sample.name or f"sample_{sample.id}"

            # Try to use bioforge_core for fast stats
            stats = {}
            try:
                import bioforge_core

                await update_job_progress(session, job, "Calculating FASTQ stats", 30)

                fq_stats = bioforge_core.calculate_stats(str(fastq_r1))
                stats["r1"] = {
                    "total_reads": fq_stats.total_reads,
                    "total_bases": fq_stats.total_bases,
                    "mean_length": fq_stats.mean_length,
                    "gc_content": fq_stats.gc_content,
                    "mean_quality": fq_stats.mean_quality,
                    "q20_percent": (
                        fq_stats.q20_bases / fq_stats.total_bases * 100
                        if fq_stats.total_bases > 0
                        else 0
                    ),
                    "q30_percent": (
                        fq_stats.q30_bases / fq_stats.total_bases * 100
                        if fq_stats.total_bases > 0
                        else 0
                    ),
                }

                if fastq_r2 and fastq_r2.exists():
                    await update_job_progress(
                        session, job, "Calculating R2 stats", 50
                    )
                    fq_stats_r2 = bioforge_core.calculate_stats(str(fastq_r2))
                    stats["r2"] = {
                        "total_reads": fq_stats_r2.total_reads,
                        "total_bases": fq_stats_r2.total_bases,
                        "mean_length": fq_stats_r2.mean_length,
                        "gc_content": fq_stats_r2.gc_content,
                        "mean_quality": fq_stats_r2.mean_quality,
                    }

            except ImportError:
                logger.warning("bioforge_core not available, using fastp")

            # Run fastp if available
            if check_tool("fastp"):
                await update_job_progress(session, job, "Running fastp QC", 60)

                html_report = output_dir / f"{sample_name}.fastp.html"
                json_report = output_dir / f"{sample_name}.fastp.json"

                fastp_cmd = [
                    "fastp",
                    "-i",
                    str(fastq_r1),
                    "-h",
                    str(html_report),
                    "-j",
                    str(json_report),
                    "--thread",
                    str(settings.pipeline.threads),
                ]

                if fastq_r2 and fastq_r2.exists():
                    fastp_cmd.extend(["-I", str(fastq_r2)])

                # Don't actually trim, just QC
                fastp_cmd.append("--disable_adapter_trimming")

                returncode, _, stderr = await run_command(fastp_cmd, timeout=1800)
                if returncode == 0 and json_report.exists():
                    import json

                    with open(json_report) as f:
                        fastp_data = json.load(f)
                        stats["fastp"] = fastp_data.get("summary", {})

            await update_job_progress(session, job, "Finalizing", 95)

            # Mark job complete
            job.status = JobStatus.COMPLETED.value
            job.progress = 100
            job.current_step = "Completed"
            job.completed_at = datetime.now(timezone.utc)
            job.result_path = str(output_dir)
            job.result_data = {"qc_stats": stats}
            await session.commit()

            logger.info(f"QC job {job_id} completed")
            return {"status": "completed", "stats": stats}

        except Exception as e:
            logger.exception(f"QC job {job_id} failed")
            job.status = JobStatus.FAILED.value
            job.error_message = str(e)
            job.completed_at = datetime.now(timezone.utc)
            await session.commit()
            return {"status": "failed", "error": str(e)}


async def run_full_pipeline(ctx: dict, job_id: int) -> dict[str, Any]:
    """
    Run full genomic analysis pipeline.

    Steps: QC -> Alignment -> Variant Calling

    Args:
        ctx: ARQ context
        job_id: ID of the master job

    Returns:
        Result dict with all outputs
    """
    logger.info(f"Starting full pipeline job {job_id}")

    results = {}

    # Step 1: QC
    qc_result = await run_qc(ctx, job_id)
    results["qc"] = qc_result
    if qc_result["status"] != "completed":
        return {"status": "failed", "step": "qc", "results": results}

    # Step 2: Alignment
    align_result = await run_alignment(ctx, job_id)
    results["alignment"] = align_result
    if align_result["status"] != "completed":
        return {"status": "failed", "step": "alignment", "results": results}

    # Step 3: Variant Calling
    vc_result = await run_variant_calling(ctx, job_id)
    results["variant_calling"] = vc_result
    if vc_result["status"] != "completed":
        return {"status": "failed", "step": "variant_calling", "results": results}

    logger.info(f"Full pipeline job {job_id} completed")
    return {"status": "completed", "results": results}


async def cleanup_old_files(ctx: dict, days_old: int = 30) -> dict[str, Any]:
    """
    Cleanup old temporary files.

    Args:
        ctx: ARQ context
        days_old: Delete files older than this many days

    Returns:
        Result dict with cleanup stats
    """
    logger.info(f"Starting cleanup of files older than {days_old} days")

    import time

    deleted_count = 0
    deleted_size = 0
    cutoff = time.time() - (days_old * 24 * 60 * 60)

    # Clean temp directory
    temp_dir = settings.storage.temp_dir
    if temp_dir.exists():
        for path in temp_dir.rglob("*"):
            if path.is_file() and path.stat().st_mtime < cutoff:
                size = path.stat().st_size
                try:
                    path.unlink()
                    deleted_count += 1
                    deleted_size += size
                except Exception as e:
                    logger.warning(f"Failed to delete {path}: {e}")

    logger.info(f"Cleanup complete: {deleted_count} files, {deleted_size / 1024 / 1024:.1f} MB")

    return {
        "status": "completed",
        "deleted_count": deleted_count,
        "deleted_size_mb": deleted_size / 1024 / 1024,
    }
