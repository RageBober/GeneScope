"""External database integrations."""

from .clinvar import ClinVarAPI
from .gnomad import GnomADAPI
from .dbsnp import DbSNPAPI
from .vcf_processor import VCFProcessor

__all__ = [
    "ClinVarAPI",
    "GnomADAPI", 
    "DbSNPAPI",
    "VCFProcessor",
]
