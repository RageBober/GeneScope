"""
PubMed Integration Module for GenoScope

Provides access to PubMed for literature search and retrieval.
"""

from __future__ import annotations

import json
import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
from xml.etree import ElementTree

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

logger = logging.getLogger(__name__)


@dataclass
class PubMedArticle:
    """Container for PubMed article information."""
    
    pmid: str
    title: str
    abstract: str
    
    # Authors
    authors: List[str] = field(default_factory=list)
    first_author: Optional[str] = None
    last_author: Optional[str] = None
    
    # Publication details
    journal: Optional[str] = None
    publication_date: Optional[str] = None
    publication_year: Optional[int] = None
    volume: Optional[str] = None
    issue: Optional[str] = None
    pages: Optional[str] = None
    
    # Identifiers
    doi: Optional[str] = None
    pmc_id: Optional[str] = None
    
    # Keywords and MeSH terms
    keywords: List[str] = field(default_factory=list)
    mesh_terms: List[str] = field(default_factory=list)
    
    # Genes and chemicals mentioned
    gene_symbols: List[str] = field(default_factory=list)
    chemicals: List[str] = field(default_factory=list)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            "pmid": self.pmid,
            "title": self.title,
            "abstract": self.abstract,
            "authors": self.authors,
            "first_author": self.first_author,
            "last_author": self.last_author,
            "journal": self.journal,
            "publication_date": self.publication_date,
            "publication_year": self.publication_year,
            "volume": self.volume,
            "issue": self.issue,
            "pages": self.pages,
            "doi": self.doi,
            "pmc_id": self.pmc_id,
            "keywords": self.keywords,
            "mesh_terms": self.mesh_terms,
            "gene_symbols": self.gene_symbols,
            "chemicals": self.chemicals
        }
    
    def get_citation(self, style: str = "apa") -> str:
        """
        Get formatted citation.
        
        Args:
            style: Citation style (apa, mla, chicago)
            
        Returns:
            Formatted citation string
        """
        if style == "apa":
            # APA style: Author, A. A. (Year). Title. Journal, Volume(Issue), pages.
            authors_str = ", ".join(self.authors[:3])
            if len(self.authors) > 3:
                authors_str += ", et al."
                
            citation = f"{authors_str} ({self.publication_year}). {self.title} "
            citation += f"{self.journal}"
            
            if self.volume:
                citation += f", {self.volume}"
            if self.issue:
                citation += f"({self.issue})"
            if self.pages:
                citation += f", {self.pages}"
                
            citation += "."
            return citation
            
        # Default to simple format
        return f"{self.authors[0] if self.authors else 'Unknown'} et al. {self.title} {self.journal}. {self.publication_year}"


class PubMedClient:
    """
    Client for accessing PubMed database.
    
    Uses NCBI E-utilities API for literature search and retrieval.
    """
    
    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self,
                 api_key: Optional[str] = None,
                 email: Optional[str] = None,
                 cache_dir: Optional[Path] = None,
                 rate_limit: float = 0.34):
        """
        Initialize PubMed client.
        
        Args:
            api_key: NCBI API key
            email: Email for NCBI (required for heavy use)
            cache_dir: Directory for caching results
            rate_limit: Minimum seconds between requests
        """
        self.api_key = api_key
        self.email = email
        self.cache_dir = Path(cache_dir) if cache_dir else Path.home() / ".genoscope" / "pubmed_cache"
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Adjust rate limit with API key
        if api_key:
            self.rate_limit = 0.1  # 10 requests/second with key
        else:
            self.rate_limit = rate_limit  # 3 requests/second without
            
        self.last_request_time = 0
        
        # Setup session
        self.session = requests.Session()
        retry_strategy = Retry(
            total=3,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504]
        )
        adapter = HTTPAdapter(max_retries=retry_strategy)
        self.session.mount("http://", adapter)
        self.session.mount("https://", adapter)
        
        # Load cache
        self._load_cache()
    
    def _load_cache(self) -> None:
        """Load cached results."""
        self.cache = {}
        cache_file = self.cache_dir / "pubmed_cache.json"
        
        if cache_file.exists():
            try:
                with open(cache_file, "r") as f:
                    self.cache = json.load(f)
                logger.info(f"Loaded {len(self.cache)} cached PubMed entries")
            except Exception as e:
                logger.warning(f"Failed to load PubMed cache: {e}")
    
    def _save_cache(self) -> None:
        """Save cache to disk."""
        cache_file = self.cache_dir / "pubmed_cache.json"
        
        try:
            with open(cache_file, "w") as f:
                json.dump(self.cache, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to save PubMed cache: {e}")
    
    def _rate_limit(self) -> None:
        """Enforce rate limiting."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        
        if time_since_last < self.rate_limit:
            time.sleep(self.rate_limit - time_since_last)
            
        self.last_request_time = time.time()
    
    def search(self,
              query: str,
              max_results: int = 20,
              sort: str = "relevance") -> List[str]:
        """
        Search PubMed for articles.
        
        Args:
            query: Search query
            max_results: Maximum number of results
            sort: Sort order (relevance, date)
            
        Returns:
            List of PMIDs
        """
        try:
            self._rate_limit()
            
            # Build search URL
            search_url = f"{self.EUTILS_BASE}/esearch.fcgi"
            params = {
                "db": "pubmed",
                "term": query,
                "retmax": max_results,
                "retmode": "json",
                "sort": sort
            }
            
            if self.api_key:
                params["api_key"] = self.api_key
            if self.email:
                params["email"] = self.email
                
            response = self.session.get(search_url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract PMIDs
            pmids = data.get("esearchresult", {}).get("idlist", [])
            
            return pmids
            
        except Exception as e:
            logger.error(f"PubMed search failed: {e}")
            return []
    
    def fetch_article(self, pmid: str) -> Optional[PubMedArticle]:
        """
        Fetch article details by PMID.
        
        Args:
            pmid: PubMed ID
            
        Returns:
            PubMedArticle object or None
        """
        # Check cache
        if pmid in self.cache:
            logger.debug(f"Using cached PubMed article {pmid}")
            return PubMedArticle(**self.cache[pmid])
        
        try:
            self._rate_limit()
            
            # Fetch article details
            fetch_url = f"{self.EUTILS_BASE}/efetch.fcgi"
            params = {
                "db": "pubmed",
                "id": pmid,
                "retmode": "xml"
            }
            
            if self.api_key:
                params["api_key"] = self.api_key
                
            response = self.session.get(fetch_url, params=params, timeout=30)
            response.raise_for_status()
            
            # Parse XML response
            article = self._parse_xml_article(response.text, pmid)
            
            if article:
                # Cache result
                self.cache[pmid] = article.to_dict()
                self._save_cache()
                
            return article
            
        except Exception as e:
            logger.error(f"Failed to fetch PubMed article {pmid}: {e}")
            return None
    
    def _parse_xml_article(self, xml_text: str, pmid: str) -> Optional[PubMedArticle]:
        """Parse PubMed XML response."""
        try:
            root = ElementTree.fromstring(xml_text)
            
            # Find article element
            article_elem = root.find(".//PubmedArticle")
            if article_elem is None:
                return None
            
            # Extract basic info
            medline = article_elem.find("MedlineCitation")
            article_data = medline.find("Article")
            
            # Title
            title = article_data.find(".//ArticleTitle")
            title_text = title.text if title is not None else ""
            
            # Abstract
            abstract_elem = article_data.find(".//Abstract/AbstractText")
            abstract_text = abstract_elem.text if abstract_elem is not None else ""
            
            # Authors
            authors = []
            author_list = article_data.find(".//AuthorList")
            if author_list is not None:
                for author in author_list.findall("Author"):
                    last_name = author.find("LastName")
                    fore_name = author.find("ForeName")
                    if last_name is not None:
                        name = last_name.text
                        if fore_name is not None:
                            name = f"{name} {fore_name.text}"
                        authors.append(name)
            
            # Journal info
            journal_elem = article_data.find(".//Journal")
            journal_name = ""
            pub_year = None
            
            if journal_elem is not None:
                journal_title = journal_elem.find("Title")
                if journal_title is not None:
                    journal_name = journal_title.text
                    
                pub_date = journal_elem.find(".//PubDate/Year")
                if pub_date is not None:
                    pub_year = int(pub_date.text)
            
            # DOI
            doi = None
            article_ids = article_elem.find(".//PubmedData/ArticleIdList")
            if article_ids is not None:
                for article_id in article_ids.findall("ArticleId"):
                    if article_id.get("IdType") == "doi":
                        doi = article_id.text
                    elif article_id.get("IdType") == "pmc":
                        pmc_id = article_id.text
            
            # MeSH terms
            mesh_terms = []
            mesh_list = medline.find(".//MeshHeadingList")
            if mesh_list is not None:
                for mesh in mesh_list.findall("MeshHeading/DescriptorName"):
                    mesh_terms.append(mesh.text)
            
            # Create article object
            article = PubMedArticle(
                pmid=pmid,
                title=title_text,
                abstract=abstract_text,
                authors=authors,
                first_author=authors[0] if authors else None,
                last_author=authors[-1] if authors else None,
                journal=journal_name,
                publication_year=pub_year,
                doi=doi,
                mesh_terms=mesh_terms
            )
            
            return article
            
        except Exception as e:
            logger.error(f"Failed to parse PubMed XML: {e}")
            return None
    
    def search_gene_variants(self,
                            gene: str,
                            variant: Optional[str] = None,
                            max_results: int = 10) -> List[PubMedArticle]:
        """
        Search for articles about gene variants.
        
        Args:
            gene: Gene symbol
            variant: Variant (optional)
            max_results: Maximum results
            
        Returns:
            List of PubMedArticle objects
        """
        # Build query
        query = f"{gene}[Gene]"
        if variant:
            query += f" AND ({variant} OR mutation OR variant)"
            
        query += " AND (pathogenic OR clinical OR disease)"
        
        # Search
        pmids = self.search(query, max_results)
        
        # Fetch articles
        articles = []
        for pmid in pmids:
            article = self.fetch_article(pmid)
            if article:
                articles.append(article)
                
        return articles
    
    def get_related_articles(self, pmid: str, max_results: int = 10) -> List[str]:
        """
        Get related articles for a PMID.
        
        Args:
            pmid: PubMed ID
            max_results: Maximum results
            
        Returns:
            List of related PMIDs
        """
        try:
            self._rate_limit()
            
            # Use elink to find related articles
            link_url = f"{self.EUTILS_BASE}/elink.fcgi"
            params = {
                "dbfrom": "pubmed",
                "db": "pubmed",
                "id": pmid,
                "cmd": "neighbor",
                "retmode": "json"
            }
            
            if self.api_key:
                params["api_key"] = self.api_key
                
            response = self.session.get(link_url, params=params, timeout=30)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract related PMIDs
            pmids = []
            link_sets = data.get("linksets", [])
            if link_sets:
                links = link_sets[0].get("linksetdbs", [])
                if links:
                    for link in links[0].get("links", [])[:max_results]:
                        pmids.append(link)
                        
            return pmids
            
        except Exception as e:
            logger.error(f"Failed to get related articles: {e}")
            return []
