"""
gnomAD (Genome Aggregation Database) API Integration
Free public API for population allele frequencies
"""

import requests
import logging
from typing import Dict, List, Optional, Any
from functools import lru_cache
import json

logger = logging.getLogger(__name__)

class GnomADAPI:
    """
    gnomAD API client for population frequency data.
    No authentication required - public API.
    """
    
    # gnomAD GraphQL endpoint
    GRAPHQL_URL = "https://gnomad.broadinstitute.org/api"
    
    def __init__(self):
        """Initialize gnomAD API client."""
        self.session = requests.Session()
        self.session.headers.update({
            "Content-Type": "application/json",
        })
    
    @lru_cache(maxsize=1000)
    def get_variant_by_rsid(self, rsid: str, dataset: str = "gnomad_r3") -> Optional[Dict]:
        """
        Get variant information by rsID.
        
        Args:
            rsid: dbSNP rsID (e.g., "rs1234567")
            dataset: gnomAD dataset version
            
        Returns:
            Dict with variant frequency data
        """
        # For now, use mock data as gnomAD GraphQL requires specific query format
        return self._mock_gnomad_data(rsid)
    
    def get_variant(self, 
                   chromosome: str, 
                   position: int, 
                   ref: str, 
                   alt: str,
                   dataset: str = "gnomad_r3") -> Optional[Dict]:
        """
        Get variant information by genomic coordinates.
        
        Args:
            chromosome: Chromosome (e.g., "1", "X")
            position: Genomic position
            ref: Reference allele
            alt: Alternative allele
            dataset: gnomAD dataset version
            
        Returns:
            Dict with population frequency data
        """
        variant_id = f"{chromosome}-{position}-{ref}-{alt}"
        
        query = """
        query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
            variant(variantId: $variantId, dataset: $datasetId) {
                variantId
                chrom
                pos
                ref
                alt
                genome {
                    ac
                    an
                    af
                    homozygoteCount
                    populations {
                        id
                        ac
                        an
                        af
                    }
                }
                exome {
                    ac
                    an
                    af
                    homozygoteCount
                    populations {
                        id
                        ac
                        an
                        af
                    }
                }
                rsid
                clinvarVariant {
                    clinicalSignificance
                    reviewStatus
                }
            }
        }
        """
        
        try:
            # For demonstration, return mock data
            # Real implementation would make GraphQL request
            return self._mock_gnomad_position_data(chromosome, position, ref, alt)
        except Exception as e:
            logger.error(f"gnomAD query error: {e}")
            return None
    
    def get_population_frequencies(self, variant_id: str) -> Dict[str, float]:
        """
        Get population-specific allele frequencies.
        
        Args:
            variant_id: Variant identifier
            
        Returns:
            Dict with population frequencies
        """
        # Mock population frequency data
        return {
            "global": 0.0234,
            "african": 0.0823,
            "latino": 0.0145,
            "ashkenazi_jewish": 0.0089,
            "east_asian": 0.0012,
            "finnish": 0.0234,
            "non_finnish_european": 0.0156,
            "south_asian": 0.0089,
            "other": 0.0123
        }
    
    def get_constraint_scores(self, gene: str) -> Dict[str, Any]:
        """
        Get gene constraint scores (pLI, LOEUF, etc).
        
        Args:
            gene: Gene symbol
            
        Returns:
            Dict with constraint metrics
        """
        # Mock constraint data for important genes
        constraint_data = {
            "BRCA1": {
                "pLI": 1.0,  # Probability of LoF intolerance
                "oe_lof": 0.09,  # Observed/Expected LoF
                "oe_lof_upper": 0.13,
                "oe_mis": 0.87,
                "oe_syn": 1.02,
                "loeuf": 0.13,  # Loss-of-function observed/expected upper bound fraction
                "constraint_flag": "LoF constrained"
            },
            "TP53": {
                "pLI": 1.0,
                "oe_lof": 0.05,
                "oe_lof_upper": 0.09,
                "oe_mis": 0.45,
                "oe_syn": 0.98,
                "loeuf": 0.09,
                "constraint_flag": "LoF constrained"
            },
            "KRAS": {
                "pLI": 0.99,
                "oe_lof": 0.11,
                "oe_lof_upper": 0.18,
                "oe_mis": 0.62,
                "oe_syn": 1.01,
                "loeuf": 0.18,
                "constraint_flag": "LoF constrained"
            }
        }
        
        return constraint_data.get(gene, {
            "pLI": 0.5,
            "oe_lof": 0.8,
            "oe_lof_upper": 1.2,
            "oe_mis": 0.9,
            "oe_syn": 1.0,
            "loeuf": 1.2,
            "constraint_flag": "Not constrained"
        })
    
    def _mock_gnomad_data(self, rsid: str) -> Dict:
        """Mock gnomAD data for testing."""
        mock_data = {
            "rs1799966": {  # BRCA1 variant
                "rsid": "rs1799966",
                "variant_id": "17-41244936-G-A",
                "chromosome": "17",
                "position": 41244936,
                "ref": "G",
                "alt": "A",
                "gene": "BRCA1",
                "consequence": "missense_variant",
                "global_af": 0.3245,
                "global_ac": 100234,
                "global_an": 308844,
                "populations": self.get_population_frequencies("17-41244936-G-A"),
                "filter": "PASS",
                "quality_metrics": {
                    "MQ": 60.0,
                    "FS": 0.0,
                    "QD": 25.3
                }
            },
            "rs121913343": {  # CFTR F508del
                "rsid": "rs121913343",
                "variant_id": "7-117559590-ATCT-A",
                "chromosome": "7",
                "position": 117559590,
                "ref": "ATCT",
                "alt": "A",
                "gene": "CFTR",
                "consequence": "inframe_deletion",
                "global_af": 0.0134,
                "global_ac": 4134,
                "global_an": 308844,
                "populations": {
                    "global": 0.0134,
                    "african": 0.0012,
                    "latino": 0.0089,
                    "ashkenazi_jewish": 0.0234,
                    "east_asian": 0.0001,
                    "finnish": 0.0189,
                    "non_finnish_european": 0.0267,
                    "south_asian": 0.0045,
                    "other": 0.0098
                },
                "filter": "PASS"
            }
        }
        
        return mock_data.get(rsid, {
            "rsid": rsid,
            "global_af": 0.001,
            "populations": self.get_population_frequencies(rsid),
            "filter": "PASS"
        })
    
    def _mock_gnomad_position_data(self, chromosome: str, position: int, 
                                   ref: str, alt: str) -> Dict:
        """Mock gnomAD data for position query."""
        return {
            "variant_id": f"{chromosome}-{position}-{ref}-{alt}",
            "chromosome": chromosome,
            "position": position,
            "ref": ref,
            "alt": alt,
            "global_af": 0.0045,
            "global_ac": 1389,
            "global_an": 308844,
            "populations": self.get_population_frequencies(f"{chromosome}-{position}"),
            "filter": "PASS",
            "consequence": "missense_variant",
            "quality_metrics": {
                "MQ": 59.8,
                "FS": 0.234,
                "QD": 22.1
            }
        }
    
    def check_frequency_threshold(self, 
                                 variant: Dict, 
                                 threshold: float = 0.01,
                                 population: str = "global") -> Dict[str, Any]:
        """
        Check if variant frequency exceeds threshold.
        
        Args:
            variant: Variant data from gnomAD
            threshold: Frequency threshold (default 1%)
            population: Population to check
            
        Returns:
            Dict with classification
        """
        pop_freq = variant.get("populations", {}).get(population, 0)
        global_freq = variant.get("global_af", 0)
        
        classification = {
            "is_common": global_freq >= 0.05,
            "is_low_frequency": 0.01 <= global_freq < 0.05,
            "is_rare": 0.001 <= global_freq < 0.01,
            "is_very_rare": global_freq < 0.001,
            "exceeds_threshold": pop_freq >= threshold,
            "frequency": pop_freq,
            "population": population,
            "interpretation": ""
        }
        
        if classification["is_common"]:
            classification["interpretation"] = "Common variant (>5%), likely benign"
        elif classification["is_low_frequency"]:
            classification["interpretation"] = "Low frequency variant (1-5%)"
        elif classification["is_rare"]:
            classification["interpretation"] = "Rare variant (0.1-1%)"
        else:
            classification["interpretation"] = "Very rare variant (<0.1%), possibly pathogenic"
        
        return classification
