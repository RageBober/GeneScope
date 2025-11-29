"""
ClinVar API Integration Module
ClinVar provides free access to information about genomic variation and its relationship to human health.
API Documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/help/
"""

import requests
import json
import time
from typing import Dict, List, Optional, Any
from datetime import datetime, timedelta
import logging
from functools import lru_cache
import xml.etree.ElementTree as ET

logger = logging.getLogger(__name__)

class ClinVarAPI:
    """
    ClinVar API client for querying variant clinical significance.
    No API key required for basic access.
    Rate limit: 3 requests per second without API key.
    """
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    CLINVAR_BASE = "https://www.ncbi.nlm.nih.gov/clinvar/variation"
    
    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        """
        Initialize ClinVar API client.
        
        Args:
            email: Email for NCBI (increases rate limit)
            api_key: NCBI API key (optional, increases rate limit to 10 req/sec)
        """
        self.email = email or "bioforge@example.com"
        self.api_key = api_key
        self.session = requests.Session()
        self.last_request_time = 0
        self.rate_limit = 0.34 if not api_key else 0.1  # seconds between requests
        
    def _rate_limit_wait(self):
        """Enforce rate limiting."""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.rate_limit:
            time.sleep(self.rate_limit - elapsed)
        self.last_request_time = time.time()
    
    def search_by_variant(self, 
                         chromosome: str, 
                         position: int, 
                         ref: str, 
                         alt: str,
                         assembly: str = "GRCh38") -> Optional[Dict]:
        """
        Search ClinVar for a specific variant.
        
        Args:
            chromosome: Chromosome (e.g., "1", "X")
            position: Genomic position
            ref: Reference allele
            alt: Alternative allele
            assembly: Genome assembly (GRCh37 or GRCh38)
            
        Returns:
            Dict with variant information or None if not found
        """
        # Build search query
        search_term = f"{chromosome}[CHR] AND {position}[CHRPOS] AND {assembly}[Assembly]"
        
        try:
            # Step 1: Search for variant IDs
            search_params = {
                "db": "clinvar",
                "term": search_term,
                "retmode": "json",
                "email": self.email
            }
            if self.api_key:
                search_params["api_key"] = self.api_key
            
            self._rate_limit_wait()
            search_response = self.session.get(
                f"{self.BASE_URL}/esearch.fcgi",
                params=search_params
            )
            search_response.raise_for_status()
            search_data = search_response.json()
            
            if not search_data.get("esearchresult", {}).get("idlist"):
                logger.info(f"No ClinVar entries found for {chromosome}:{position} {ref}>{alt}")
                return None
            
            # Step 2: Fetch variant details
            variant_ids = search_data["esearchresult"]["idlist"]
            return self.get_variant_details(variant_ids[0])
            
        except requests.RequestException as e:
            logger.error(f"ClinVar API error: {e}")
            return None
    
    def get_variant_details(self, variant_id: str) -> Optional[Dict]:
        """
        Get detailed information about a ClinVar variant.
        
        Args:
            variant_id: ClinVar variation ID
            
        Returns:
            Dict with variant details
        """
        try:
            # Fetch variant summary
            summary_params = {
                "db": "clinvar",
                "id": variant_id,
                "retmode": "xml",
                "email": self.email
            }
            if self.api_key:
                summary_params["api_key"] = self.api_key
            
            self._rate_limit_wait()
            response = self.session.get(
                f"{self.BASE_URL}/esummary.fcgi",
                params=summary_params
            )
            response.raise_for_status()
            
            # Parse XML response
            root = ET.fromstring(response.text)
            
            variant_info = {
                "clinvar_id": variant_id,
                "url": f"{self.CLINVAR_BASE}/{variant_id}/",
                "clinical_significance": None,
                "review_status": None,
                "conditions": [],
                "gene": None,
                "last_evaluated": None,
                "variation_type": None,
                "molecular_consequence": None
            }
            
            # Extract information from XML
            doc_sum = root.find(".//DocumentSummary")
            if doc_sum is not None:
                # Clinical significance
                clin_sig = doc_sum.find(".//clinical_significance")
                if clin_sig is not None:
                    variant_info["clinical_significance"] = clin_sig.text
                
                # Gene symbol
                gene = doc_sum.find(".//gene_symbol")
                if gene is not None:
                    variant_info["gene"] = gene.text
                
                # Conditions
                trait_set = doc_sum.find(".//trait_set")
                if trait_set is not None:
                    for trait in trait_set.findall(".//trait_name"):
                        if trait.text:
                            variant_info["conditions"].append(trait.text)
                
                # Review status
                review = doc_sum.find(".//review_status")
                if review is not None:
                    variant_info["review_status"] = review.text
                
                # Variation type
                var_type = doc_sum.find(".//variation_type")
                if var_type is not None:
                    variant_info["variation_type"] = var_type.text
            
            return variant_info
            
        except Exception as e:
            logger.error(f"Error fetching variant details: {e}")
            return None
    
    def get_clinical_significance_summary(self, variant_id: str) -> str:
        """
        Get a human-readable summary of clinical significance.
        
        Args:
            variant_id: ClinVar variation ID
            
        Returns:
            String summary of clinical significance
        """
        details = self.get_variant_details(variant_id)
        if not details:
            return "No clinical data available"
        
        significance = details.get("clinical_significance", "Unknown")
        conditions = details.get("conditions", [])
        review = details.get("review_status", "Not reviewed")
        
        summary = f"Clinical Significance: {significance}\n"
        if conditions:
            summary += f"Associated Conditions: {', '.join(conditions[:3])}\n"
        summary += f"Review Status: {review}"
        
        return summary
    
    @lru_cache(maxsize=1000)
    def search_by_rsid(self, rsid: str) -> Optional[Dict]:
        """
        Search ClinVar by dbSNP rsID.
        
        Args:
            rsid: dbSNP rsID (e.g., "rs1234567")
            
        Returns:
            Dict with variant information
        """
        if not rsid.startswith("rs"):
            rsid = f"rs{rsid}"
        
        search_term = f"{rsid}[RS]"
        
        try:
            search_params = {
                "db": "clinvar",
                "term": search_term,
                "retmode": "json",
                "email": self.email
            }
            if self.api_key:
                search_params["api_key"] = self.api_key
            
            self._rate_limit_wait()
            response = self.session.get(
                f"{self.BASE_URL}/esearch.fcgi",
                params=search_params
            )
            response.raise_for_status()
            data = response.json()
            
            if data.get("esearchresult", {}).get("idlist"):
                return self.get_variant_details(data["esearchresult"]["idlist"][0])
            
            return None
            
        except Exception as e:
            logger.error(f"Error searching by rsID {rsid}: {e}")
            return None
    
    def batch_query_variants(self, variants: List[Dict]) -> List[Dict]:
        """
        Query multiple variants in batch.
        
        Args:
            variants: List of dicts with keys: chromosome, position, ref, alt
            
        Returns:
            List of variant information dicts
        """
        results = []
        for variant in variants:
            result = self.search_by_variant(
                chromosome=variant["chromosome"],
                position=variant["position"],
                ref=variant["ref"],
                alt=variant["alt"],
                assembly=variant.get("assembly", "GRCh38")
            )
            if result:
                result["query"] = variant
            results.append(result)
            
        return results
    
    def interpret_significance(self, significance: str) -> Dict[str, Any]:
        """
        Interpret ClinVar clinical significance into categories.
        
        Args:
            significance: ClinVar significance string
            
        Returns:
            Dict with interpretation details
        """
        pathogenic_terms = ["pathogenic", "likely pathogenic"]
        benign_terms = ["benign", "likely benign"]
        uncertain_terms = ["uncertain significance", "vus", "vous"]
        
        sig_lower = significance.lower() if significance else ""
        
        interpretation = {
            "category": "unknown",
            "is_pathogenic": False,
            "is_benign": False,
            "is_uncertain": False,
            "requires_follow_up": False,
            "actionable": False
        }
        
        if any(term in sig_lower for term in pathogenic_terms):
            interpretation.update({
                "category": "pathogenic",
                "is_pathogenic": True,
                "requires_follow_up": True,
                "actionable": True
            })
        elif any(term in sig_lower for term in benign_terms):
            interpretation.update({
                "category": "benign",
                "is_benign": True,
                "requires_follow_up": False,
                "actionable": False
            })
        elif any(term in sig_lower for term in uncertain_terms):
            interpretation.update({
                "category": "uncertain",
                "is_uncertain": True,
                "requires_follow_up": True,
                "actionable": False
            })
        elif "conflicting" in sig_lower:
            interpretation.update({
                "category": "conflicting",
                "requires_follow_up": True,
                "actionable": False
            })
        
        return interpretation


# Example usage and testing
if __name__ == "__main__":
    # Initialize ClinVar API
    clinvar = ClinVarAPI()
    
    # Test 1: Search for a known pathogenic variant (BRCA1)
    print("Testing ClinVar API...")
    print("-" * 50)
    
    # Example: BRCA1 pathogenic variant
    result = clinvar.search_by_variant(
        chromosome="17",
        position=43044295,
        ref="G",
        alt="A",
        assembly="GRCh38"
    )
    
    if result:
        print(f"Found variant: {result['clinvar_id']}")
        print(f"Clinical Significance: {result['clinical_significance']}")
        print(f"Gene: {result['gene']}")
        print(f"Conditions: {', '.join(result['conditions'][:3])}")
        print(f"ClinVar URL: {result['url']}")
    else:
        print("No result found")
    
    print("\n" + "-" * 50)
    
    # Test 2: Search by rsID
    rs_result = clinvar.search_by_rsid("rs121913343")
    if rs_result:
        print(f"Found by rsID: {rs_result['clinical_significance']}")
    
    print("\n" + "-" * 50)
    
    # Test 3: Interpret significance
    interpretation = clinvar.interpret_significance("Pathogenic")
    print(f"Interpretation: {interpretation}")
