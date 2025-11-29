# src/genoscope/pipeline/real_genomics_pipeline.py

import os
import sys
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any
import logging
from dataclasses import dataclass
import json
import hashlib
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import shutil
from datetime import datetime

import pysam
from cyvcf2 import VCF, Writer
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)

@dataclass
class PipelineConfig:
    """Конфигурация pipeline"""
    # Пути к инструментам
    bwa_path: str = "bwa"
    samtools_path: str = "samtools"
    bcftools_path: str = "bcftools"
    gatk_path: str = "gatk"
    fastp_path: str = "fastp"
    fastqc_path: str = "fastqc"
    snpeff_path: str = "snpEff"
    vep_path: str = "vep"
    
    # Параметры обработки
    threads: int = 4
    memory_gb: int = 8
    quality_threshold: int = 20
    min_coverage: int = 10
    
    # Пути к референсам
    reference_genome: str = "/data/references/GRCh38/GRCh38.fa"
    known_sites_vcf: str = "/data/references/dbsnp_156.vcf.gz"
    cosmic_vcf: str = "/data/references/cosmic_v99.vcf.gz"
    clinvar_vcf: str = "/data/references/clinvar_20240101.vcf.gz"

class RealGenomicsPipeline:
    """Производственный pipeline для обработки геномных данных"""
    
    def __init__(self, config: Optional[PipelineConfig] = None):
        self.config = config or PipelineConfig()
        self.work_dir = Path(tempfile.mkdtemp(prefix="genoscope_"))
        self._validate_tools()
        self._index_reference()
        
    def _validate_tools(self) -> bool:
        """Проверка доступности инструментов"""
        tools_status = {}
        
        critical_tools = [
            self.config.bwa_path,
            self.config.samtools_path,
            self.config.bcftools_path,
            self.config.gatk_path
        ]
        
        for tool in critical_tools:
            tools_status[tool] = shutil.which(tool) is not None
            if not tools_status[tool]:
                logger.warning(f"Tool {tool} not found in PATH")
        
        # Проверяем наличие референсного генома
        if not Path(self.config.reference_genome).exists():
            logger.error(f"Reference genome not found: {self.config.reference_genome}")
            return False
            
        return all(tools_status.values())
    
    def _index_reference(self):
        """Индексация референсного генома если нужно"""
        ref_path = Path(self.config.reference_genome)
        
        # BWA индекс
        if not (ref_path.parent / f"{ref_path.name}.bwt").exists():
            logger.info("Creating BWA index...")
            subprocess.run([
                self.config.bwa_path, "index", 
                str(ref_path)
            ])
        
        # Samtools индекс
        if not (ref_path.parent / f"{ref_path.name}.fai").exists():
            logger.info("Creating samtools index...")
            subprocess.run([
                self.config.samtools_path, "faidx",
                str(ref_path)
            ])
        
        # GATK sequence dictionary
        dict_file = ref_path.parent / f"{ref_path.stem}.dict"
        if not dict_file.exists():
            logger.info("Creating sequence dictionary...")
            subprocess.run([
                self.config.gatk_path, "CreateSequenceDictionary",
                f"-R={ref_path}",
                f"-O={dict_file}"
            ])
    
    def run_full_pipeline(self,
                         sample_name: str,
                         fastq_r1: str,
                         fastq_r2: Optional[str] = None,
                         output_dir: str = "./results") -> Dict[str, Any]:
        """
        Запуск полного pipeline от FASTQ до аннотированного VCF
        
        Этапы:
        1. QC и препроцессинг
        2. Выравнивание
        3. Обработка BAM
        4. Variant calling
        5. Фильтрация
        6. Аннотация
        """
        
        results = {
            "sample": sample_name,
            "start_time": datetime.now().isoformat(),
            "stages": {}
        }
        
        output_path = Path(output_dir) / sample_name
        output_path.mkdir(parents=True, exist_ok=True)
        
        try:
            # 1. QC и препроцессинг
            logger.info("Stage 1: QC and preprocessing")
            clean_r1, clean_r2 = self._run_qc_and_trimming(
                fastq_r1, fastq_r2, output_path / "qc"
            )
            results["stages"]["qc"] = "completed"
            
            # 2. Выравнивание
            logger.info("Stage 2: Alignment")
            aligned_bam = self._run_alignment(
                clean_r1, clean_r2, sample_name, output_path / "alignment"
            )
            results["stages"]["alignment"] = "completed"
            
            # 3. Обработка BAM
            logger.info("Stage 3: BAM processing")
            processed_bam = self._process_bam(
                aligned_bam, sample_name, output_path / "processed"
            )
            results["stages"]["bam_processing"] = "completed"
            
            # 4. Variant calling
            logger.info("Stage 4: Variant calling")
            raw_vcf = self._call_variants(
                processed_bam, sample_name, output_path / "variants"
            )
            results["stages"]["variant_calling"] = "completed"
            
            # 5. Фильтрация вариантов
            logger.info("Stage 5: Variant filtering")
            filtered_vcf = self._filter_variants(
                raw_vcf, output_path / "filtered"
            )
            results["stages"]["filtering"] = "completed"
            
            # 6. Аннотация
            logger.info("Stage 6: Annotation")
            annotated_vcf = self._annotate_variants(
                filtered_vcf, output_path / "annotated"
            )
            results["stages"]["annotation"] = "completed"
            
            # 7. Генерация отчета
            logger.info("Stage 7: Report generation")
            report = self._generate_report(
                annotated_vcf, output_path / "report"
            )
            results["stages"]["report"] = "completed"
            
            results["output"] = {
                "final_vcf": str(annotated_vcf),
                "report": str(report),
                "bam": str(processed_bam)
            }
            
        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            results["error"] = str(e)
            results["status"] = "failed"
        else:
            results["status"] = "completed"
        finally:
            results["end_time"] = datetime.now().isoformat()
            
        return results
    
    def _run_qc_and_trimming(self, 
                            fastq_r1: str,
                            fastq_r2: Optional[str],
                            output_dir: Path) -> Tuple[str, Optional[str]]:
        """QC и очистка ридов"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # FastQC до обработки
        subprocess.run([
            self.config.fastqc_path,
            "-o", str(output_dir),
            "-t", str(self.config.threads),
            fastq_r1
        ])
        
        # Fastp для очистки
        clean_r1 = output_dir / "clean_R1.fastq.gz"
        clean_r2 = output_dir / "clean_R2.fastq.gz" if fastq_r2 else None
        
        fastp_cmd = [
            self.config.fastp_path,
            "-i", fastq_r1,
            "-o", str(clean_r1),
            "-h", str(output_dir / "fastp.html"),
            "-j", str(output_dir / "fastp.json"),
            "--thread", str(self.config.threads),
            "--qualified_quality_phred", str(self.config.quality_threshold),
            "--length_required", "50",
            "--detect_adapter_for_pe"
        ]
        
        if fastq_r2:
            fastp_cmd.extend(["-I", fastq_r2, "-O", str(clean_r2)])
        
        subprocess.run(fastp_cmd, check=True)
        
        return str(clean_r1), str(clean_r2) if clean_r2 else None
    
    def _run_alignment(self,
                      fastq_r1: str,
                      fastq_r2: Optional[str],
                      sample_name: str,
                      output_dir: Path) -> str:
        """Выравнивание на референсный геном"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        aligned_sam = output_dir / f"{sample_name}.sam"
        aligned_bam = output_dir / f"{sample_name}.bam"
        sorted_bam = output_dir / f"{sample_name}.sorted.bam"
        
        # BWA MEM выравнивание
        bwa_cmd = [
            self.config.bwa_path, "mem",
            "-t", str(self.config.threads),
            "-R", f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA",
            self.config.reference_genome,
            fastq_r1
        ]
        
        if fastq_r2:
            bwa_cmd.append(fastq_r2)
        
        # Выравнивание и конвертация в BAM
        with open(aligned_sam, 'w') as sam_file:
            subprocess.run(bwa_cmd, stdout=sam_file, check=True)
        
        # SAM to BAM
        subprocess.run([
            self.config.samtools_path, "view",
            "-b", "-@", str(self.config.threads),
            "-o", str(aligned_bam),
            str(aligned_sam)
        ], check=True)
        
        # Сортировка
        subprocess.run([
            self.config.samtools_path, "sort",
            "-@", str(self.config.threads),
            "-o", str(sorted_bam),
            str(aligned_bam)
        ], check=True)
        
        # Индексация
        subprocess.run([
            self.config.samtools_path, "index",
            str(sorted_bam)
        ], check=True)
        
        # Удаляем промежуточные файлы
        os.remove(aligned_sam)
        os.remove(aligned_bam)
        
        return str(sorted_bam)
    
    def _process_bam(self,
                    bam_file: str,
                    sample_name: str,
                    output_dir: Path) -> str:
        """Обработка BAM файла - mark duplicates, BQSR"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        marked_bam = output_dir / f"{sample_name}.marked.bam"
        recal_bam = output_dir / f"{sample_name}.recal.bam"
        
        # Mark duplicates с GATK
        subprocess.run([
            self.config.gatk_path, "MarkDuplicates",
            f"-I={bam_file}",
            f"-O={marked_bam}",
            f"-M={output_dir}/metrics.txt",
            "--CREATE_INDEX=true"
        ], check=True)
        
        # Base Quality Score Recalibration (BQSR)
        recal_table = output_dir / "recal_data.table"
        
        # Создаем таблицу рекалибровки
        subprocess.run([
            self.config.gatk_path, "BaseRecalibrator",
            f"-R={self.config.reference_genome}",
            f"-I={marked_bam}",
            f"--known-sites={self.config.known_sites_vcf}",
            f"-O={recal_table}"
        ], check=True)
        
        # Применяем рекалибровку
        subprocess.run([
            self.config.gatk_path, "ApplyBQSR",
            f"-R={self.config.reference_genome}",
            f"-I={marked_bam}",
            f"--bqsr-recal-file={recal_table}",
            f"-O={recal_bam}"
        ], check=True)
        
        return str(recal_bam)
    
    def _call_variants(self,
                      bam_file: str,
                      sample_name: str,
                      output_dir: Path) -> str:
        """Вызов вариантов с помощью GATK HaplotypeCaller"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        raw_vcf = output_dir / f"{sample_name}.raw.vcf.gz"
        
        subprocess.run([
            self.config.gatk_path, "HaplotypeCaller",
            f"-R={self.config.reference_genome}",
            f"-I={bam_file}",
            f"-O={raw_vcf}",
            "--emit-ref-confidence", "GVCF",
            "--native-pair-hmm-threads", str(self.config.threads)
        ], check=True)
        
        return str(raw_vcf)
    
    def _filter_variants(self,
                        vcf_file: str,
                        output_dir: Path) -> str:
        """Фильтрация вариантов по качеству"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        filtered_vcf = output_dir / Path(vcf_file).name.replace(".raw.", ".filtered.")
        
        # Фильтрация с bcftools
        subprocess.run([
            self.config.bcftools_path, "filter",
            "-i", f"QUAL>{self.config.quality_threshold} && DP>{self.config.min_coverage}",
            "-o", str(filtered_vcf),
            "-O", "z",  # compressed output
            vcf_file
        ], check=True)
        
        # Индексация
        subprocess.run([
            self.config.bcftools_path, "index",
            str(filtered_vcf)
        ], check=True)
        
        return str(filtered_vcf)
    
    def _annotate_variants(self,
                          vcf_file: str,
                          output_dir: Path) -> str:
        """Аннотация вариантов с SnpEff/VEP"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        annotated_vcf = output_dir / Path(vcf_file).name.replace(".filtered.", ".annotated.")
        
        # SnpEff аннотация
        subprocess.run([
            "java", "-jar", self.config.snpeff_path,
            "GRCh38.99",  # версия базы данных
            vcf_file,
            "-v",
            "-stats", str(output_dir / "snpeff_summary.html"),
            ">", str(annotated_vcf)
        ], shell=True, check=True)
        
        return str(annotated_vcf)
    
    def _generate_report(self,
                        vcf_file: str,
                        output_dir: Path) -> Path:
        """Генерация HTML отчета"""
        output_dir.mkdir(parents=True, exist_ok=True)
        
        report_path = output_dir / "genomic_report.html"
        
        # Анализ VCF файла
        stats = self._analyze_vcf(vcf_file)
        
        # Генерация HTML
        html_content = self._create_html_report(stats)
        
        report_path.write_text(html_content)
        
        return report_path
    
    def _analyze_vcf(self, vcf_file: str) -> Dict[str, Any]:
        """Анализ VCF для отчета"""
        vcf = VCF(vcf_file)
        
        stats = {
            "total_variants": 0,
            "snps": 0,
            "indels": 0,
            "transitions": 0,
            "transversions": 0,
            "heterozygous": 0,
            "homozygous": 0,
            "quality_distribution": [],
            "depth_distribution": []
        }
        
        for variant in vcf:
            stats["total_variants"] += 1
            
            # Тип варианта
            if len(variant.REF) == 1 and len(variant.ALT[0]) == 1:
                stats["snps"] += 1
                
                # Transitions vs Transversions
                if (variant.REF in "AG" and variant.ALT[0] in "AG") or \
                   (variant.REF in "CT" and variant.ALT[0] in "CT"):
                    stats["transitions"] += 1
                else:
                    stats["transversions"] += 1
            else:
                stats["indels"] += 1
            
            # Генотип
            gt = variant.gt_types[0] if variant.gt_types else None
            if gt == 1:  # 0/1
                stats["heterozygous"] += 1
            elif gt == 3:  # 1/1
                stats["homozygous"] += 1
            
            # Качество и глубина
            stats["quality_distribution"].append(variant.QUAL)
            if variant.INFO.get("DP"):
                stats["depth_distribution"].append(variant.INFO["DP"])
        
        vcf.close()
        
        # Вычисляем Ti/Tv ratio
        if stats["transversions"] > 0:
            stats["titv_ratio"] = stats["transitions"] / stats["transversions"]
        else:
            stats["titv_ratio"] = 0
        
        return stats
    
    def _create_html_report(self, stats: Dict[str, Any]) -> str:
        """Создание HTML отчета"""
        html = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>Genomic Analysis Report</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 40px; }}
                h1 {{ color: #2c3e50; }}
                .stat-card {{ 
                    background: #f8f9fa; 
                    border-radius: 8px; 
                    padding: 20px; 
                    margin: 20px 0;
                }}
                .metric {{ 
                    display: inline-block; 
                    margin: 10px 20px; 
                }}
                .metric-value {{ 
                    font-size: 24px; 
                    font-weight: bold; 
                    color: #3498db; 
                }}
                .metric-label {{ 
                    color: #7f8c8d; 
                    font-size: 14px; 
                }}
            </style>
        </head>
        <body>
            <h1>🧬 Genomic Analysis Report</h1>
            
            <div class="stat-card">
                <h2>Variant Statistics</h2>
                <div class="metric">
                    <div class="metric-value">{stats['total_variants']:,}</div>
                    <div class="metric-label">Total Variants</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{stats['snps']:,}</div>
                    <div class="metric-label">SNPs</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{stats['indels']:,}</div>
                    <div class="metric-label">InDels</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{stats.get('titv_ratio', 0):.2f}</div>
                    <div class="metric-label">Ti/Tv Ratio</div>
                </div>
            </div>
            
            <div class="stat-card">
                <h2>Genotype Distribution</h2>
                <div class="metric">
                    <div class="metric-value">{stats['heterozygous']:,}</div>
                    <div class="metric-label">Heterozygous</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{stats['homozygous']:,}</div>
                    <div class="metric-label">Homozygous</div>
                </div>
            </div>
            
            <div class="stat-card">
                <h2>Quality Metrics</h2>
                <div class="metric">
                    <div class="metric-value">
                        {np.mean(stats['quality_distribution']):.1f}
                    </div>
                    <div class="metric-label">Mean Quality</div>
                </div>
                <div class="metric">
                    <div class="metric-value">
                        {np.mean(stats['depth_distribution']) if stats['depth_distribution'] else 0:.1f}
                    </div>
                    <div class="metric-label">Mean Depth</div>
                </div>
            </div>
            
            <div style="margin-top: 40px; color: #95a5a6; font-size: 12px;">
                Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </div>
        </body>
        </html>
        """
        
        return html


# Использование
if __name__ == "__main__":
    # Создаем конфигурацию
    config = PipelineConfig(
        threads=8,
        memory_gb=16,
        reference_genome="/data/references/GRCh38/GRCh38.fa"
    )
    
    # Инициализируем pipeline
    pipeline = RealGenomicsPipeline(config)
    
    # Запускаем обработку
    results = pipeline.run_full_pipeline(
        sample_name="PATIENT_001",
        fastq_r1="/data/samples/PATIENT_001_R1.fastq.gz",
        fastq_r2="/data/samples/PATIENT_001_R2.fastq.gz",
        output_dir="/data/results"
    )
    
    print(json.dumps(results, indent=2))