"""
🔬 Taxonomy Analysis Utilities

Provides helper functions for taxonomic data processing:
- Taxonomic lineage parsing
- Taxa filtering and normalization
- Cross-tool taxonomy mapping
"""

from __future__ import annotations

import logging
from collections import defaultdict
from typing import Any, Optional

logger = logging.getLogger(__name__)


class TaxonomyAnalyzer:
    """
    Utility class for taxonomic analysis and normalization.

    Examples:
        >>> analyzer = TaxonomyAnalyzer()
        >>> lineage = analyzer.parse_lineage("Bacteria;Proteobacteria;Escherichia coli")
        >>> print(lineage)
        ['Bacteria', 'Proteobacteria', 'Escherichia coli']
    """

    # Standard taxonomic ranks (in order)
    RANKS = [
        "domain", "kingdom", "phylum", "class", "order",
        "family", "genus", "species", "strain"
    ]

    # Rank abbreviations used by different tools
    RANK_ABBREV = {
        "D": "domain",
        "K": "kingdom",
        "P": "phylum",
        "C": "class",
        "O": "order",
        "F": "family",
        "G": "genus",
        "S": "species",
        "T": "strain",
    }

    def __init__(self):
        """Initialize taxonomy analyzer"""
        pass

    def parse_lineage(
        self,
        lineage_str: str,
        separator: str = ";",
    ) -> list[str]:
        """
        Parse taxonomic lineage string.

        Args:
            lineage_str: Lineage string (e.g., "Bacteria;Proteobacteria;E. coli")
            separator: Separator character

        Returns:
            List of taxonomic names

        Examples:
            >>> analyzer = TaxonomyAnalyzer()
            >>> analyzer.parse_lineage("Bacteria;Firmicutes;Bacillus")
            ['Bacteria', 'Firmicutes', 'Bacillus']
        """
        if not lineage_str:
            return []

        # Split and clean
        taxa = [t.strip() for t in lineage_str.split(separator)]

        # Remove empty entries
        taxa = [t for t in taxa if t]

        return taxa

    def normalize_rank_name(self, rank_abbrev: str) -> str:
        """
        Normalize rank abbreviation to full name.

        Args:
            rank_abbrev: Rank abbreviation (e.g., 'S', 'G')

        Returns:
            Full rank name (e.g., 'species', 'genus')

        Examples:
            >>> analyzer = TaxonomyAnalyzer()
            >>> analyzer.normalize_rank_name('S')
            'species'
            >>> analyzer.normalize_rank_name('genus')
            'genus'
        """
        # If already full name, return as-is
        if rank_abbrev.lower() in self.RANKS:
            return rank_abbrev.lower()

        # Convert abbreviation
        return self.RANK_ABBREV.get(rank_abbrev.upper(), rank_abbrev.lower())

    def filter_by_rank(
        self,
        taxa_list: list[dict[str, Any]],
        rank: str,
    ) -> list[dict[str, Any]]:
        """
        Filter taxa by taxonomic rank.

        Args:
            taxa_list: List of taxa dictionaries (must have 'rank' key)
            rank: Target rank (e.g., 'species', 'genus', or 'S', 'G')

        Returns:
            Filtered list of taxa

        Examples:
            >>> taxa = [
            ...     {'rank': 'S', 'name': 'E. coli'},
            ...     {'rank': 'G', 'name': 'Escherichia'},
            ... ]
            >>> analyzer = TaxonomyAnalyzer()
            >>> species = analyzer.filter_by_rank(taxa, 'S')
            >>> len(species)
            1
        """
        normalized_rank = self.normalize_rank_name(rank)

        # Also check for abbreviations
        rank_variants = {rank, normalized_rank}
        for abbrev, full_name in self.RANK_ABBREV.items():
            if full_name == normalized_rank:
                rank_variants.add(abbrev)

        return [
            t for t in taxa_list
            if t.get("rank", "").lower() in rank_variants
            or t.get("rank", "") in rank_variants
        ]

    def aggregate_by_rank(
        self,
        taxa_list: list[dict[str, Any]],
        rank: str,
    ) -> dict[str, int]:
        """
        Aggregate taxa counts by rank.

        Args:
            taxa_list: List of taxa dictionaries
            rank: Target rank to aggregate

        Returns:
            Dictionary mapping taxon name to total count

        Examples:
            >>> taxa = [
            ...     {'rank': 'S', 'name': 'E. coli', 'count': 100},
            ...     {'rank': 'S', 'name': 'E. coli', 'count': 50},
            ...     {'rank': 'S', 'name': 'S. aureus', 'count': 75},
            ... ]
            >>> analyzer = TaxonomyAnalyzer()
            >>> counts = analyzer.aggregate_by_rank(taxa, 'S')
            >>> counts['E. coli']
            150
        """
        filtered = self.filter_by_rank(taxa_list, rank)

        aggregated = defaultdict(int)

        for taxon in filtered:
            name = taxon.get("name", "Unknown")
            count = taxon.get("num_taxon_reads", taxon.get("count", 1))
            aggregated[name] += count

        return dict(aggregated)

    def get_top_taxa(
        self,
        taxa_list: list[dict[str, Any]],
        n: int = 10,
        rank: Optional[str] = None,
        sort_by: str = "percentage",
    ) -> list[dict[str, Any]]:
        """
        Get top N taxa from classification results.

        Args:
            taxa_list: List of taxa dictionaries
            n: Number of top taxa to return
            rank: Filter by rank (optional)
            sort_by: Sort field ('percentage', 'count', etc.)

        Returns:
            Top N taxa

        Examples:
            >>> taxa = [
            ...     {'name': 'E. coli', 'percentage': 45.0},
            ...     {'name': 'S. aureus', 'percentage': 30.0},
            ... ]
            >>> analyzer = TaxonomyAnalyzer()
            >>> top = analyzer.get_top_taxa(taxa, n=1)
            >>> top[0]['name']
            'E. coli'
        """
        # Filter by rank if specified
        if rank:
            taxa_list = self.filter_by_rank(taxa_list, rank)

        # Sort by specified field
        if sort_by in taxa_list[0] if taxa_list else False:
            sorted_taxa = sorted(
                taxa_list,
                key=lambda x: x.get(sort_by, 0),
                reverse=True,
            )
        else:
            sorted_taxa = taxa_list

        return sorted_taxa[:n]


# ═══════════════════════════════════════════════════════════════
#  Convenience functions
# ═══════════════════════════════════════════════════════════════

def parse_lineage(lineage_str: str, separator: str = ";") -> list[str]:
    """
    Parse taxonomic lineage string (convenience function).

    Args:
        lineage_str: Lineage string
        separator: Separator character

    Returns:
        List of taxonomic names
    """
    analyzer = TaxonomyAnalyzer()
    return analyzer.parse_lineage(lineage_str, separator)


def normalize_rank(rank: str) -> str:
    """
    Normalize rank abbreviation (convenience function).

    Args:
        rank: Rank abbreviation or full name

    Returns:
        Full rank name
    """
    analyzer = TaxonomyAnalyzer()
    return analyzer.normalize_rank_name(rank)
