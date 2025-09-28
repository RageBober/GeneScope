"""
Structural Variant Analysis Module
Handles CNVs, translocations, inversions, and other structural variants
"""

import logging
from typing import Dict, List, Optional, Any, Tuple
from dataclasses import dataclass, field
from enum import Enum
import numpy as np

logger = logging.getLogger(__name__)

class SVType(Enum):
    """Types of structural variants."""
    DELETION = "DEL"
    DUPLICATION = "DUP"
    INVERSION = "INV"
    TRANSLOCATION = "TRA"
    INSERTION = "INS"
    CNV_GAIN = "CNV+"
    CNV_LOSS = "CNV-"
    COMPLEX = "CPX"

@dataclass
class StructuralVariant:
    """Represents a structural variant."""
    sv_type: SVType
    chromosome: str
    start: int
    end: int
    size: int
    chromosome2: Optional[str] = None  # For translocations
    start2: Optional[int] = None
    end2: Optional[int] = None
    copy_number: Optional[float] = None
    confidence: float = 0.0
    supporting_reads: int = 0
    genes_affected: List[str] = field(default_factory=list)
    
    def __post_init__(self):
        # Calculate size if not provided
        if self.size == 0 and self.end and self.start:
            self.size = abs(self.end - self.start)
    
    def is_pathogenic_size(self) -> bool:
        """Check if SV size suggests pathogenicity."""
        if self.sv_type in [SVType.DELETION, SVType.DUPLICATION]:
            return self.size > 1000000  # >1Mb
        return False
    
    def get_classification(self) -> str:
        """Classify SV by size."""
        if self.size < 50:
            return "Small (<50bp)"
        elif self.size < 1000:
            return "Medium (50bp-1kb)"
        elif self.size < 1000000:
            return "Large (1kb-1Mb)"
        else:
            return "Very Large (>1Mb)"


class StructuralVariantAnalyzer:
    """Analyzes structural variants for clinical significance."""
    
    def __init__(self):
        """Initialize SV analyzer."""
        self.pathogenic_regions = self._load_pathogenic_regions()
        self.dosage_sensitive_genes = self._load_dosage_sensitive_genes()
    
    def _load_pathogenic_regions(self) -> List[Dict]:
        """Load known pathogenic CNV regions."""
        # Examples of well-known pathogenic CNV regions
        return [
            {
                "name": "DiGeorge syndrome",
                "chr": "22",
                "start": 18900000,
                "end": 21500000,
                "type": "deletion",
                "phenotype": "Heart defects, immune deficiency, developmental delay"
            },
            {
                "name": "Williams-Beuren syndrome",
                "chr": "7",
                "start": 72700000,
                "end": 74100000,
                "type": "deletion",
                "phenotype": "Cardiovascular disease, developmental delay, distinctive facies"
            },
            {
                "name": "Prader-Willi/Angelman",
                "chr": "15",
                "start": 22700000,
                "end": 28100000,
                "type": "deletion",
                "phenotype": "Hypotonia, intellectual disability, behavioral problems"
            },
            {
                "name": "1p36 deletion",
                "chr": "1",
                "start": 1,
                "end": 5400000,
                "type": "deletion",
                "phenotype": "Intellectual disability, seizures, hearing loss"
            }
        ]
    
    def _load_dosage_sensitive_genes(self) -> Dict[str, Dict]:
        """Load genes sensitive to copy number changes."""
        return {
            "SHANK3": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "Phelan-McDermid"},
            "NRXN1": {"pHaplo": 0.95, "pTriplo": 0.05, "syndrome": "Autism, schizophrenia"},
            "CNTNAP2": {"pHaplo": 0.90, "pTriplo": 0.10, "syndrome": "Autism, epilepsy"},
            "CHD7": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "CHARGE syndrome"},
            "TBX1": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "DiGeorge syndrome"},
            "ELN": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "Williams-Beuren"},
            "MECP2": {"pHaplo": 0.99, "pTriplo": 0.99, "syndrome": "Rett syndrome"},
            "FMR1": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "Fragile X"},
            "PTEN": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "PTEN hamartoma"},
            "NF1": {"pHaplo": 0.99, "pTriplo": 0.01, "syndrome": "Neurofibromatosis"}
        }
    
    def analyze_cnv(self, cnv: StructuralVariant) -> Dict[str, Any]:
        """
        Analyze a copy number variant for pathogenicity.
        
        Args:
            cnv: CNV to analyze
            
        Returns:
            Analysis results with interpretation
        """
        results = {
            "variant": cnv,
            "classification": "Uncertain",
            "evidence": [],
            "affected_regions": [],
            "dosage_sensitive_genes": [],
            "size_category": cnv.get_classification(),
            "interpretation": ""
        }
        
        # Check size
        if cnv.is_pathogenic_size():
            results["evidence"].append(f"Large CNV ({cnv.size/1000000:.1f}Mb)")
            results["classification"] = "Likely Pathogenic"
        
        # Check overlap with known pathogenic regions
        for region in self.pathogenic_regions:
            if self._overlaps_region(cnv, region):
                results["affected_regions"].append(region["name"])
                results["evidence"].append(f"Overlaps {region['name']} region")
                results["classification"] = "Pathogenic"
                results["interpretation"] = region["phenotype"]
        
        # Check dosage-sensitive genes
        for gene in cnv.genes_affected:
            if gene in self.dosage_sensitive_genes:
                gene_info = self.dosage_sensitive_genes[gene]
                results["dosage_sensitive_genes"].append(gene)
                
                if cnv.sv_type in [SVType.DELETION, SVType.CNV_LOSS]:
                    if gene_info["pHaplo"] > 0.9:
                        results["evidence"].append(f"{gene} is haploinsufficient")
                        results["classification"] = "Pathogenic"
                elif cnv.sv_type in [SVType.DUPLICATION, SVType.CNV_GAIN]:
                    if gene_info["pTriplo"] > 0.9:
                        results["evidence"].append(f"{gene} is triplosensitive")
                        results["classification"] = "Pathogenic"
        
        # Generate interpretation
        if not results["interpretation"]:
            results["interpretation"] = self._generate_interpretation(results)
        
        return results
    
    def _overlaps_region(self, cnv: StructuralVariant, region: Dict) -> bool:
        """Check if CNV overlaps with a genomic region."""
        if cnv.chromosome != region["chr"]:
            return False
        
        # Check for overlap
        return not (cnv.end < region["start"] or cnv.start > region["end"])
    
    def _generate_interpretation(self, analysis: Dict) -> str:
        """Generate clinical interpretation for CNV."""
        classification = analysis["classification"]
        
        if classification == "Pathogenic":
            interp = "This CNV is classified as pathogenic. "
            if analysis["affected_regions"]:
                interp += f"It affects known syndrome regions: {', '.join(analysis['affected_regions'])}. "
            if analysis["dosage_sensitive_genes"]:
                interp += f"Dosage-sensitive genes affected: {', '.join(analysis['dosage_sensitive_genes'])}. "
            interp += "Clinical correlation and genetic counseling are recommended."
        
        elif classification == "Likely Pathogenic":
            interp = "This CNV is likely pathogenic based on size and gene content. "
            interp += "Further clinical evaluation is recommended."
        
        else:
            interp = "The clinical significance of this CNV is uncertain. "
            interp += "Consider family studies and clinical correlation."
        
        return interp
    
    def detect_complex_rearrangements(self, variants: List[StructuralVariant]) -> List[Dict]:
        """
        Detect complex genomic rearrangements from multiple SVs.
        
        Args:
            variants: List of structural variants
            
        Returns:
            List of detected complex events
        """
        complex_events = []
        
        # Chromothripsis detection (multiple rearrangements in single chromosome)
        chr_counts = {}
        for sv in variants:
            chr_counts[sv.chromosome] = chr_counts.get(sv.chromosome, 0) + 1
        
        for chr, count in chr_counts.items():
            if count >= 10:
                complex_events.append({
                    "type": "Chromothripsis",
                    "chromosome": chr,
                    "sv_count": count,
                    "interpretation": f"Possible chromothripsis on chromosome {chr}"
                })
        
        # Chromoplexy detection (interconnected translocations)
        translocations = [sv for sv in variants if sv.sv_type == SVType.TRANSLOCATION]
        if len(translocations) >= 3:
            connected = self._find_connected_translocations(translocations)
            if connected:
                complex_events.append({
                    "type": "Chromoplexy",
                    "chromosomes": list(set([t.chromosome for t in connected] + 
                                          [t.chromosome2 for t in connected if t.chromosome2])),
                    "sv_count": len(connected),
                    "interpretation": "Complex interconnected translocations detected"
                })
        
        return complex_events
    
    def _find_connected_translocations(self, translocations: List[StructuralVariant]) -> List[StructuralVariant]:
        """Find interconnected translocations."""
        # Simplified: return all if multiple chromosomes involved
        chromosomes = set()
        for t in translocations:
            chromosomes.add(t.chromosome)
            if t.chromosome2:
                chromosomes.add(t.chromosome2)
        
        if len(chromosomes) >= 3:
            return translocations
        return []
    
    def calculate_tumor_purity(self, cnvs: List[StructuralVariant]) -> float:
        """
        Estimate tumor purity from CNV data.
        
        Args:
            cnvs: List of CNVs with copy number data
            
        Returns:
            Estimated tumor purity (0-1)
        """
        # Simplified calculation based on deviation from diploid
        if not cnvs:
            return 0.0
        
        deviations = []
        for cnv in cnvs:
            if cnv.copy_number is not None:
                # Deviation from diploid (CN=2)
                dev = abs(cnv.copy_number - 2.0)
                deviations.append(dev)
        
        if deviations:
            # Estimate purity from maximum deviation
            max_dev = max(deviations)
            # Assuming pure tumor would have integer copy numbers
            purity = min(max_dev / 2.0, 1.0)
            return purity
        
        return 0.0


def parse_sv_from_vcf_info(info_field: str) -> Optional[StructuralVariant]:
    """
    Parse structural variant from VCF INFO field.
    
    Args:
        info_field: VCF INFO field string
        
    Returns:
        StructuralVariant object or None
    """
    # Parse INFO field for SV information
    info_dict = {}
    for item in info_field.split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            info_dict[key] = value
    
    # Determine SV type
    svtype = info_dict.get('SVTYPE', '').upper()
    if not svtype:
        return None
    
    try:
        sv_type_map = {
            'DEL': SVType.DELETION,
            'DUP': SVType.DUPLICATION,
            'INV': SVType.INVERSION,
            'TRA': SVType.TRANSLOCATION,
            'INS': SVType.INSERTION,
            'CNV': SVType.CNV_GAIN
        }
        
        sv_type = sv_type_map.get(svtype, SVType.COMPLEX)
        
        # Parse coordinates
        chr1 = info_dict.get('CHR', '1')
        start = int(info_dict.get('POS', 0))
        end = int(info_dict.get('END', start + 1))
        size = int(info_dict.get('SVLEN', end - start))
        
        # For translocations
        chr2 = info_dict.get('CHR2')
        
        # Supporting evidence
        supporting_reads = int(info_dict.get('SR', 0)) + int(info_dict.get('PE', 0))
        
        return StructuralVariant(
            sv_type=sv_type,
            chromosome=chr1,
            start=start,
            end=end,
            size=abs(size),
            chromosome2=chr2,
            supporting_reads=supporting_reads,
            confidence=float(info_dict.get('QUAL', 0))
        )
        
    except (ValueError, KeyError) as e:
        logger.error(f"Error parsing SV: {e}")
        return None
