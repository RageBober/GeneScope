"""
dbSNP API Integration Module
Integrates with NCBI's dbSNP database for SNP information.
"""

import requests
import logging
import time
from typing import Dict, List, Optional, Any
from functools import lru_cache
import json

logger = logging.getLogger(__name__)

class DbSNPAPI:
    """
    dbSNP API client for querying SNP information.
    Uses NCBI E-utilities API (free, no key required).
    """
    
    BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self, email: Optional[str] = None, api_key: Optional[str] = None):
        """
        Initialize dbSNP API client.
        
        Args:
            email: Email for NCBI (required for E-utilities)
            api_key: NCBI API key (optional, increases rate limit)
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
    
    @lru_cache(maxsize=1000)
    def get_snp_info(self, rsid: str) -> Optional[Dict]:
        """
        Get information about a SNP by rsID.
        
        Args:
            rsid: dbSNP rsID (e.g., "rs1234567" or just "1234567")
            
        Returns:
            Dict with SNP information or None if not found
        """
        # Clean rsID
        if not rsid.startswith("rs"):
            rsid = f"rs{rsid}"
        
        # Remove 'rs' for the API query
        snp_id = rsid[2:]
        
        try:
            # Search for SNP
            search_params = {
                "db": "snp",
                "term": rsid,
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
                logger.info(f"No dbSNP entry found for {rsid}")
                return None
            
            # Get SNP details
            snp_uid = search_data["esearchresult"]["idlist"][0]
            
            summary_params = {
                "db": "snp",
                "id": snp_uid,
                "retmode": "json",
                "email": self.email
            }
            if self.api_key:
                summary_params["api_key"] = self.api_key
            
            self._rate_limit_wait()
            summary_response = self.session.get(
                f"{self.BASE_URL}/esummary.fcgi",
                params=summary_params
            )
            summary_response.raise_for_status()
            summary_data = summary_response.json()
            
            # Parse SNP information
            if "result" in summary_data and snp_uid in summary_data["result"]:
                snp_data = summary_data["result"][snp_uid]
                
                return {
                    "rsid": rsid,
                    "uid": snp_uid,
                    "chromosome": snp_data.get("chr", ""),
                    "position": snp_data.get("chrpos", ""),
                    "gene": snp_data.get("gene_name", ""),
                    "alleles": snp_data.get("allele", ""),
                    "maf": snp_data.get("global_maf", None),  # Minor allele frequency
                    "clinical_significance": snp_data.get("clinical_significance", ""),
                    "validated": snp_data.get("validated", False),
                    "url": f"https://www.ncbi.nlm.nih.gov/snp/{rsid}"
                }
            
            return None
            
        except Exception as e:
            logger.error(f"Error fetching SNP info for {rsid}: {e}")
            return None
    
    def batch_query_snps(self, rsids: List[str]) -> List[Optional[Dict]]:
        """
        Query multiple SNPs in batch.
        
        Args:
            rsids: List of rsIDs
            
        Returns:
            List of SNP information dicts (None for not found)
        """
        results = []
        for rsid in rsids:
            result = self.get_snp_info(rsid)
            results.append(result)
        return results
    
    def get_population_frequencies(self, rsid: str) -> Optional[Dict]:
        """
        Get population-specific allele frequencies for a SNP.
        
        Args:
            rsid: dbSNP rsID
            
        Returns:
            Dict with population frequencies
        """
        snp_info = self.get_snp_info(rsid)
        if not snp_info:
            return None
        
        # This would typically require additional API calls to get detailed population data
        # For now, return basic MAF if available
        return {
            "rsid": rsid,
            "global_maf": snp_info.get("maf"),
            "populations": {
                "global": snp_info.get("maf", "N/A")
            }
        }
