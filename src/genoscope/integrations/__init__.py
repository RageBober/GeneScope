"""
GenoScope Integrations Module
External database and service integrations for genomic analysis.
"""

from .clinvar import ClinVarAPI
from .dbsnp import DbSNPAPI
from .vcf_processor import VCFProcessor, Variant
from .cosmic import COSMICAPI
from .gnomad import GnomADAPI

__all__ = [
    "ClinVarAPI",
    "DbSNPAPI", 
    "VCFProcessor",
    "Variant",
    "COSMICAPI",
    "GnomADAPI"
]
