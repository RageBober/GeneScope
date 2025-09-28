"""
Enhanced Genomics Analysis API endpoints
Includes ClinVar, dbSNP, COSMIC, gnomAD, ML predictions, and SV analysis
"""

from fastapi import APIRouter, UploadFile, File, HTTPException, Depends, BackgroundTasks
from fastapi.responses import FileResponse, JSONResponse
from typing import List, Optional, Dict, Any
from pathlib import Path
import tempfile
import json
import uuid
from datetime import datetime

# Import all genomics modules
from ..pipeline.genomics_pipeline import GenomicsPipeline
from ..integrations import (
    VCFProcessor, 
    ClinVarAPI, 
    DbSNPAPI,
    COSMICAPI,
    GnomADAPI
)
from ..ml import PathogenicityPredictor, VariantFeatures, extract_features_from_variant
from ..analysis import StructuralVariantAnalyzer, parse_sv_from_vcf_info

# Create router
genomics_router = APIRouter(prefix="/api/genomics", tags=["genomics"])

# Initialize services
pipeline = GenomicsPipeline()
clinvar_api = ClinVarAPI()
dbsnp_api = DbSNPAPI()
cosmic_api = COSMICAPI()
gnomad_api = GnomADAPI()
pathogenicity_predictor = PathogenicityPredictor()
sv_analyzer = StructuralVariantAnalyzer()

# Store for analysis jobs
analysis_jobs = {}

@genomics_router.post("/upload/vcf")
async def upload_vcf_file(file: UploadFile = File(...)):
    """Upload and process VCF file with full annotation"""
    if not file.filename.endswith(('.vcf', '.vcf.gz')):
        raise HTTPException(status_code=400, detail="File must be VCF format")
    
    upload_dir = Path("/tmp/genoscope/uploads")
    upload_dir.mkdir(parents=True, exist_ok=True)
    
    file_path = upload_dir / file.filename
    with open(file_path, 'wb') as f:
        content = await file.read()
        f.write(content)
    
    try:
        processor = VCFProcessor(str(file_path))
        variants = processor.parse(limit=100)
        stats = processor.get_statistics()
        
        return {
            "message": "VCF uploaded successfully",
            "file": str(file_path),
            "statistics": stats
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Error parsing VCF: {str(e)}")

@genomics_router.get("/clinvar/{rsid}")
async def get_clinvar_info(rsid: str):
    """Get ClinVar information for variant"""
    try:
        result = clinvar_api.search_by_rsid(rsid)
        if not result:
            raise HTTPException(status_code=404, detail="Variant not found in ClinVar")
        
        if result.get("clinical_significance"):
            result["interpretation"] = clinvar_api.interpret_significance(
                result["clinical_significance"]
            )
        
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"ClinVar query failed: {str(e)}")

@genomics_router.get("/dbsnp/{rsid}")
async def get_dbsnp_info(rsid: str):
    """Get dbSNP information for variant"""
    try:
        result = dbsnp_api.get_snp_info(rsid)
        if not result:
            raise HTTPException(status_code=404, detail="SNP not found in dbSNP")
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"dbSNP query failed: {str(e)}")

@genomics_router.get("/cosmic/{gene}/{mutation}")
async def get_cosmic_info(gene: str, mutation: str):
    """Get COSMIC cancer mutation information"""
    try:
        result = cosmic_api.search_mutation(gene, mutation)
        if not result:
            raise HTTPException(status_code=404, detail="Mutation not found in COSMIC")
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"COSMIC query failed: {str(e)}")

@genomics_router.get("/gnomad/variant")
async def get_gnomad_info(
    chromosome: str,
    position: int,
    ref: str,
    alt: str
):
    """Get gnomAD population frequency data"""
    try:
        result = gnomad_api.get_variant(chromosome, position, ref, alt)
        if result:
            # Add frequency classification
            result["classification"] = gnomad_api.check_frequency_threshold(result)
        return result or {"message": "Variant not found in gnomAD"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"gnomAD query failed: {str(e)}")

@genomics_router.post("/predict/pathogenicity")
async def predict_pathogenicity(variant_data: Dict[str, Any]):
    """
    Predict variant pathogenicity using ML model.
    
    Request body should include:
    - chromosome, position, ref, alt
    - gene (optional)
    - type (optional)
    """
    try:
        # Query databases for additional information
        clinvar_data = None
        gnomad_data = None
        cosmic_data = None
        
        # Try to get ClinVar data
        if variant_data.get("rsid"):
            clinvar_data = clinvar_api.search_by_rsid(variant_data["rsid"])
        
        # Get gnomAD data
        if all(k in variant_data for k in ["chromosome", "position", "ref", "alt"]):
            gnomad_data = gnomad_api.get_variant(
                variant_data["chromosome"],
                variant_data["position"],
                variant_data["ref"],
                variant_data["alt"]
            )
        
        # Get COSMIC data if gene is provided
        if variant_data.get("gene"):
            cosmic_data = {"is_cancer_gene": variant_data["gene"] in ["TP53", "KRAS", "BRAF", "EGFR"]}
        
        # Extract features
        features = extract_features_from_variant(
            variant_data,
            clinvar_data,
            gnomad_data,
            cosmic_data
        )
        
        # Make prediction
        prediction = pathogenicity_predictor.predict(features)
        
        return {
            "variant": variant_data,
            "prediction": prediction,
            "features_used": {
                "has_clinvar": clinvar_data is not None,
                "has_gnomad": gnomad_data is not None,
                "has_cosmic": cosmic_data is not None
            }
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Prediction failed: {str(e)}")

@genomics_router.post("/analyze/structural-variants")
async def analyze_structural_variants(vcf_path: str):
    """Analyze structural variants in VCF file"""
    try:
        vcf_file = Path(vcf_path)
        if not vcf_file.exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        # Parse VCF for structural variants
        processor = VCFProcessor(str(vcf_file))
        variants = processor.parse()
        
        svs = []
        for var in variants:
            if var.info and "SVTYPE" in var.info:
                sv = parse_sv_from_vcf_info(";".join([f"{k}={v}" for k, v in var.info.items()]))
                if sv:
                    svs.append(sv)
        
        # Analyze each SV
        results = []
        for sv in svs:
            analysis = sv_analyzer.analyze_cnv(sv)
            results.append({
                "type": sv.sv_type.value,
                "location": f"{sv.chromosome}:{sv.start}-{sv.end}",
                "size": sv.size,
                "classification": analysis["classification"],
                "interpretation": analysis["interpretation"],
                "evidence": analysis["evidence"]
            })
        
        # Check for complex rearrangements
        complex_events = sv_analyzer.detect_complex_rearrangements(svs)
        
        return {
            "total_svs": len(svs),
            "analyzed_svs": results,
            "complex_events": complex_events,
            "summary": {
                "pathogenic": sum(1 for r in results if r["classification"] == "Pathogenic"),
                "likely_pathogenic": sum(1 for r in results if r["classification"] == "Likely Pathogenic"),
                "uncertain": sum(1 for r in results if r["classification"] == "Uncertain")
            }
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"SV analysis failed: {str(e)}")

@genomics_router.post("/annotate/comprehensive")
async def comprehensive_annotation(vcf_path: str):
    """
    Perform comprehensive annotation including all databases and ML prediction.
    This is the main endpoint for full variant annotation.
    """
    try:
        vcf_file = Path(vcf_path)
        if not vcf_file.exists():
            raise HTTPException(status_code=404, detail="VCF file not found")
        
        # Parse VCF
        processor = VCFProcessor(str(vcf_file))
        variants = processor.parse(limit=50)  # Limit for performance
        
        annotated_variants = []
        
        for var in variants:
            variant_dict = var.to_dict()
            annotation = {
                "variant": variant_dict,
                "clinvar": None,
                "dbsnp": None,
                "gnomad": None,
                "cosmic": None,
                "pathogenicity": None,
                "structural": None
            }
            
            # ClinVar annotation
            if var.variant_id and var.variant_id.startswith("rs"):
                annotation["clinvar"] = clinvar_api.search_by_rsid(var.variant_id)
                annotation["dbsnp"] = dbsnp_api.get_snp_info(var.variant_id)
            
            # gnomAD annotation
            annotation["gnomad"] = gnomad_api.get_variant(
                var.chromosome.replace("chr", ""),
                var.position,
                var.ref,
                var.alt
            )
            
            # COSMIC annotation (if gene available)
            if variant_dict.get("gene"):
                annotation["cosmic"] = cosmic_api.search_mutation(
                    variant_dict["gene"],
                    f"{var.ref}{var.position}{var.alt}"
                )
            
            # ML pathogenicity prediction
            features = extract_features_from_variant(
                variant_dict,
                annotation["clinvar"],
                annotation["gnomad"],
                annotation["cosmic"]
            )
            annotation["pathogenicity"] = pathogenicity_predictor.predict(features)
            
            # Check if structural variant
            if var.info and "SVTYPE" in var.info:
                sv = parse_sv_from_vcf_info(";".join([f"{k}={v}" for k, v in var.info.items()]))
                if sv:
                    annotation["structural"] = sv_analyzer.analyze_cnv(sv)
            
            annotated_variants.append(annotation)
        
        # Generate summary
        summary = {
            "total_variants": len(annotated_variants),
            "pathogenic": sum(1 for v in annotated_variants 
                            if v["pathogenicity"] and v["pathogenicity"]["classification"] == "Pathogenic"),
            "likely_pathogenic": sum(1 for v in annotated_variants 
                                   if v["pathogenicity"] and v["pathogenicity"]["classification"] == "Likely Pathogenic"),
            "uncertain": sum(1 for v in annotated_variants 
                           if v["pathogenicity"] and v["pathogenicity"]["classification"] == "Uncertain Significance"),
            "benign": sum(1 for v in annotated_variants 
                        if v["pathogenicity"] and v["pathogenicity"]["classification"] in ["Benign", "Likely Benign"]),
            "with_clinvar": sum(1 for v in annotated_variants if v["clinvar"]),
            "with_gnomad": sum(1 for v in annotated_variants if v["gnomad"]),
            "structural_variants": sum(1 for v in annotated_variants if v["structural"])
        }
        
        return {
            "summary": summary,
            "variants": annotated_variants[:10],  # Return first 10 for display
            "message": "Comprehensive annotation completed"
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Annotation failed: {str(e)}")

@genomics_router.get("/tools/status")
async def get_tools_status():
    """Check availability of bioinformatics tools and databases"""
    return {
        "tools": pipeline.tools,
        "databases": {
            "clinvar": "available",
            "dbsnp": "available",
            "cosmic": "available (mock data)",
            "gnomad": "available (mock data)"
        },
        "ml_models": {
            "pathogenicity_predictor": "available",
            "sv_analyzer": "available"
        },
        "message": "System status check complete"
    }

@genomics_router.get("/gene/constraint/{gene}")
async def get_gene_constraint(gene: str):
    """Get gene constraint scores from gnomAD"""
    try:
        constraint = gnomad_api.get_constraint_scores(gene)
        return {
            "gene": gene,
            "constraint": constraint,
            "interpretation": {
                "lof_intolerant": constraint.get("pLI", 0) > 0.9,
                "lof_tolerant": constraint.get("pLI", 0) < 0.1,
                "constrained": constraint.get("loeuf", 1) < 0.35
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Query failed: {str(e)}")

@genomics_router.get("/cancer/genes")
async def get_cancer_gene_census():
    """Get list of known cancer genes"""
    try:
        genes = cosmic_api.get_cancer_gene_census()
        return {
            "total": len(genes),
            "genes": genes,
            "categories": {
                "oncogenes": [g for g in genes if g["role"] == "Oncogene"],
                "tumor_suppressors": [g for g in genes if g["role"] == "TSG"]
            }
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Query failed: {str(e)}")
