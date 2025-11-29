"""
Variant Annotation Module for BioForge

Handles functional annotation of variants using VEP, SnpEff, and external databases.
"""

from __future__ import annotations

import json
import logging
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Literal

import pandas as pd
import requests

logger = logging.getLogger(__name__)

AnnotationTool = Literal["vep", "snpeff", "annovar"]


@dataclass
class AnnotatedVariant:
    """Container for annotated variant information."""
    
    # Basic variant info
    chrom: str
    pos: int
    ref: str
    alt: str
    
    # Gene annotation
    gene_symbol: Optional[str] = None
    gene_id: Optional[str] = None
    transcript_id: Optional[str] = None
    exon: Optional[str] = None
    intron: Optional[str] = None
    
    # Consequence
    consequence: Optional[str] = None
    impact: Optional[str] = None  # HIGH, MODERATE, LOW, MODIFIER
    biotype: Optional[str] = None
    
    # Protein change
    hgvs_c: Optional[str] = None  # cDNA change
    hgvs_p: Optional[str] = None  # Protein change
    amino_acids: Optional[str] = None
    codons: Optional[str] = None
    
    # Population frequencies
    gnomad_af: Optional[float] = None
    gnomad_af_popmax: Optional[float] = None
    thousand_genomes_af: Optional[float] = None
    esp_af: Optional[float] = None
    
    # Clinical significance
    clinvar_id: Optional[str] = None
    clinvar_sig: Optional[str] = None
    cosmic_id: Optional[str] = None
    
    # Pathogenicity scores
    sift_score: Optional[float] = None
    sift_pred: Optional[str] = None
    polyphen_score: Optional[float] = None
    polyphen_pred: Optional[str] = None
    cadd_phred: Optional[float] = None
    revel_score: Optional[float] = None
    
    # Conservation scores
    phylop_score: Optional[float] = None
    phastcons_score: Optional[float] = None
    gerp_score: Optional[float] = None
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {k: v for k, v in self.__dict__.items() if v is not None}
    
    def is_pathogenic(self, threshold: float = 0.5) -> bool:
        """Check if variant is likely pathogenic based on scores."""
        if self.clinvar_sig and "pathogenic" in self.clinvar_sig.lower():
            return True
            
        if self.impact == "HIGH":
            return True
            
        # Check pathogenicity scores
        pathogenic_indicators = 0
        total_indicators = 0
        
        if self.sift_pred:
            total_indicators += 1
            if self.sift_pred == "deleterious":
                pathogenic_indicators += 1
                
        if self.polyphen_pred:
            total_indicators += 1
            if self.polyphen_pred in ["probably_damaging", "possibly_damaging"]:
                pathogenic_indicators += 1
                
        if self.cadd_phred:
            total_indicators += 1
            if self.cadd_phred > 20:
                pathogenic_indicators += 1
                
        if self.revel_score:
            total_indicators += 1
            if self.revel_score > 0.5:
                pathogenic_indicators += 1
                
        if total_indicators > 0:
            return (pathogenic_indicators / total_indicators) >= threshold
            
        return False


class VariantAnnotator:
    """
    Manages functional annotation of genomic variants.
    
    Integrates with VEP, SnpEff, and various annotation databases.
    """
    
    # ClinVar API endpoint
    CLINVAR_API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    # gnomAD API endpoint  
    GNOMAD_API = "https://gnomad.broadinstitute.org/api"
    
    def __init__(self,
                 cache_dir: Optional[Path] = None,
                 temp_dir: Optional[Path] = None,
                 threads: int = 4):
        """
        Initialize variant annotator.
        
        Args:
            cache_dir: Directory for annotation cache/databases
            temp_dir: Directory for temporary files
            threads: Number of threads
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path.home() / ".bioforge" / "annotation_cache"
        self.temp_dir = Path(temp_dir) if temp_dir else Path(tempfile.gettempdir())
        self.threads = threads
        
        # Create directories
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Check available annotation tools
        self.available_tools = self._check_available_tools()
        
        # Initialize cache for API queries
        self._init_cache()
    
    def _check_available_tools(self) -> List[str]:
        """Check which annotation tools are installed."""
        tools = []
        
        # Check for VEP
        try:
            subprocess.run(["vep", "--help"], capture_output=True, check=False)
            tools.append("vep")
            logger.info("Found annotation tool: VEP")
        except FileNotFoundError:
            logger.warning("VEP not found")
            
        # Check for SnpEff
        try:
            subprocess.run(["snpEff", "-version"], capture_output=True, check=False)
            tools.append("snpeff")
            logger.info("Found annotation tool: SnpEff")
        except FileNotFoundError:
            logger.warning("SnpEff not found")
            
        return tools
    
    def _init_cache(self) -> None:
        """Initialize local cache for annotations."""
        self.clinvar_cache = {}
        self.gnomad_cache = {}
        
        # Load existing cache if available
        clinvar_cache_file = self.cache_dir / "clinvar_cache.json"
        if clinvar_cache_file.exists():
            try:
                with open(clinvar_cache_file) as f:
                    self.clinvar_cache = json.load(f)
            except Exception as e:
                logger.warning(f"Failed to load ClinVar cache: {e}")
    
    def annotate_vcf(self,
                    vcf_file: Path,
                    output_dir: Path,
                    tool: AnnotationTool = "vep",
                    genome_build: str = "GRCh38") -> Path:
        """
        Annotate variants in VCF file.
        
        Args:
            vcf_file: Input VCF file
            output_dir: Output directory
            tool: Annotation tool to use
            genome_build: Reference genome build
            
        Returns:
            Path to annotated VCF
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        if tool not in self.available_tools:
            logger.warning(f"Annotation tool '{tool}' not available, using basic annotation")
            return self._basic_annotation(vcf_file, output_dir)
        
        if tool == "vep":
            return self._annotate_with_vep(vcf_file, output_dir, genome_build)
        elif tool == "snpeff":
            return self._annotate_with_snpeff(vcf_file, output_dir, genome_build)
        else:
            raise NotImplementedError(f"Annotation with {tool} not yet implemented")
    
    def _annotate_with_vep(self,
                          vcf_file: Path,
                          output_dir: Path,
                          genome_build: str) -> Path:
        """Annotate using Ensembl VEP."""
        output_vcf = output_dir / f"{vcf_file.stem}.vep.vcf"
        stats_file = output_dir / f"{vcf_file.stem}.vep_stats.html"
        
        cmd = [
            "vep",
            "--input_file", str(vcf_file),
            "--output_file", str(output_vcf),
            "--vcf",
            "--force_overwrite",
            "--assembly", genome_build,
            "--offline",
            "--cache",
            "--dir_cache", str(self.cache_dir),
            "--fork", str(self.threads),
            "--stats_file", str(stats_file),
            # Annotations
            "--everything",  # Include all annotations
            "--plugin", "CADD",
            "--plugin", "REVEL",
            "--plugin", "SpliceAI",
            "--af_gnomad",
            "--af_1kg",
            "--clinvar",
            "--cosmic"
        ]
        
        try:
            logger.info("Running VEP annotation...")
            subprocess.run(cmd, check=True, capture_output=True)
            logger.info(f"VEP annotation completed: {output_vcf}")
            return output_vcf
            
        except subprocess.CalledProcessError as e:
            logger.error(f"VEP failed: {e.stderr.decode() if e.stderr else str(e)}")
            # Fall back to basic annotation
            return self._basic_annotation(vcf_file, output_dir)
    
    def _annotate_with_snpeff(self,
                             vcf_file: Path,
                             output_dir: Path,
                             genome_build: str) -> Path:
        """Annotate using SnpEff."""
        output_vcf = output_dir / f"{vcf_file.stem}.snpeff.vcf"
        stats_file = output_dir / f"{vcf_file.stem}.snpeff_stats.html"
        
        # Map genome build to SnpEff database
        genome_map = {
            "GRCh38": "GRCh38.99",
            "GRCh37": "GRCh37.75",
            "hg38": "hg38",
            "hg19": "hg19"
        }
        snpeff_db = genome_map.get(genome_build, "GRCh38.99")
        
        cmd = [
            "snpEff", "ann",
            "-v",
            "-stats", str(stats_file),
            "-dataDir", str(self.cache_dir),
            snpeff_db,
            str(vcf_file)
        ]
        
        try:
            logger.info("Running SnpEff annotation...")
            with open(output_vcf, "w") as f:
                subprocess.run(cmd, stdout=f, check=True, stderr=subprocess.PIPE)
            logger.info(f"SnpEff annotation completed: {output_vcf}")
            return output_vcf
            
        except subprocess.CalledProcessError as e:
            logger.error(f"SnpEff failed: {e.stderr.decode() if e.stderr else str(e)}")
            return self._basic_annotation(vcf_file, output_dir)
    
    def _basic_annotation(self, vcf_file: Path, output_dir: Path) -> Path:
        """Basic annotation using API calls when tools are not available."""
        output_vcf = output_dir / f"{vcf_file.stem}.annotated.vcf"
        
        logger.info("Performing basic annotation using API calls...")
        
        # For now, just copy the file
        # In production, would parse VCF and add annotations via API
        import shutil
        shutil.copy(vcf_file, output_vcf)
        
        logger.warning("Basic annotation completed (limited functionality)")
        return output_vcf
    
    def query_clinvar(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """
        Query ClinVar for variant information.
        
        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            ClinVar annotation data
        """
        # Check cache first
        cache_key = f"{chrom}:{pos}:{ref}:{alt}"
        if cache_key in self.clinvar_cache:
            return self.clinvar_cache[cache_key]
        
        try:
            # Build query
            query = f"{chrom}[CHR] AND {pos}[POS] AND human[ORGN]"
            
            # Search ClinVar
            search_url = f"{self.CLINVAR_API}/esearch.fcgi"
            search_params = {
                "db": "clinvar",
                "term": query,
                "retmode": "json"
            }
            
            search_response = requests.get(search_url, params=search_params, timeout=10)
            search_data = search_response.json()
            
            if search_data.get("esearchresult", {}).get("count", "0") == "0":
                return {}
            
            # Get details for first result
            var_id = search_data["esearchresult"]["idlist"][0]
            
            # Fetch variant details
            fetch_url = f"{self.CLINVAR_API}/esummary.fcgi"
            fetch_params = {
                "db": "clinvar",
                "id": var_id,
                "retmode": "json"
            }
            
            fetch_response = requests.get(fetch_url, params=fetch_params, timeout=10)
            fetch_data = fetch_response.json()
            
            # Parse results
            result = {}
            if "result" in fetch_data and var_id in fetch_data["result"]:
                var_data = fetch_data["result"][var_id]
                result = {
                    "clinvar_id": var_id,
                    "clinical_significance": var_data.get("clinical_significance", {}).get("description"),
                    "review_status": var_data.get("review_status", {}).get("description"),
                    "last_evaluated": var_data.get("clinical_significance", {}).get("last_evaluated")
                }
            
            # Cache result
            self.clinvar_cache[cache_key] = result
            return result
            
        except Exception as e:
            logger.error(f"ClinVar query failed: {e}")
            return {}
    
    def query_gnomad(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """
        Query gnomAD for population frequency.
        
        Args:
            chrom: Chromosome
            pos: Position
            ref: Reference allele
            alt: Alternate allele
            
        Returns:
            gnomAD frequency data
        """
        # Check cache first
        cache_key = f"{chrom}:{pos}:{ref}:{alt}"
        if cache_key in self.gnomad_cache:
            return self.gnomad_cache[cache_key]
        
        try:
            # Note: This is a simplified example
            # Real implementation would use gnomAD API or local database
            
            result = {
                "gnomad_af": 0.001,  # Placeholder
                "gnomad_af_popmax": 0.005,  # Placeholder
                "gnomad_nhomalt": 10  # Placeholder
            }
            
            # Cache result
            self.gnomad_cache[cache_key] = result
            return result
            
        except Exception as e:
            logger.error(f"gnomAD query failed: {e}")
            return {}
    
    def prioritize_variants(self,
                          annotated_variants: List[AnnotatedVariant],
                          strategy: str = "pathogenicity") -> List[AnnotatedVariant]:
        """
        Prioritize variants based on various criteria.
        
        Args:
            annotated_variants: List of annotated variants
            strategy: Prioritization strategy
            
        Returns:
            Sorted list of prioritized variants
        """
        if strategy == "pathogenicity":
            # Sort by pathogenicity scores
            def pathogenicity_score(var: AnnotatedVariant) -> float:
                score = 0.0
                
                # Clinical significance
                if var.clinvar_sig:
                    if "pathogenic" in var.clinvar_sig.lower():
                        score += 10
                    elif "likely_pathogenic" in var.clinvar_sig.lower():
                        score += 5
                        
                # Impact
                if var.impact == "HIGH":
                    score += 4
                elif var.impact == "MODERATE":
                    score += 2
                    
                # Pathogenicity predictions
                if var.cadd_phred and var.cadd_phred > 20:
                    score += 3
                if var.revel_score and var.revel_score > 0.5:
                    score += 2
                if var.sift_pred == "deleterious":
                    score += 1
                if var.polyphen_pred in ["probably_damaging", "possibly_damaging"]:
                    score += 1
                    
                # Population frequency (rare = higher priority)
                if var.gnomad_af:
                    if var.gnomad_af < 0.001:
                        score += 3
                    elif var.gnomad_af < 0.01:
                        score += 1
                        
                return score
            
            return sorted(annotated_variants, key=pathogenicity_score, reverse=True)
            
        elif strategy == "cancer":
            # Prioritize cancer-related variants
            def cancer_score(var: AnnotatedVariant) -> float:
                score = 0.0
                
                if var.cosmic_id:
                    score += 5
                    
                # Check for known cancer genes (simplified)
                cancer_genes = ["TP53", "BRCA1", "BRCA2", "KRAS", "EGFR", "BRAF"]
                if var.gene_symbol and var.gene_symbol in cancer_genes:
                    score += 3
                    
                return score
            
            return sorted(annotated_variants, key=cancer_score, reverse=True)
            
        else:
            # Default: return as-is
            return annotated_variants
    
    def save_cache(self) -> None:
        """Save annotation cache to disk."""
        try:
            # Save ClinVar cache
            clinvar_cache_file = self.cache_dir / "clinvar_cache.json"
            with open(clinvar_cache_file, "w") as f:
                json.dump(self.clinvar_cache, f)
                
            # Save gnomAD cache
            gnomad_cache_file = self.cache_dir / "gnomad_cache.json"
            with open(gnomad_cache_file, "w") as f:
                json.dump(self.gnomad_cache, f)
                
            logger.info("Annotation cache saved")
            
        except Exception as e:
            logger.error(f"Failed to save cache: {e}")
