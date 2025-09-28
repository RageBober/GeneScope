"""
Machine Learning Module for Genomics Analysis
"""

from .pathogenicity_model import (
    PathogenicityPredictor,
    VariantFeatures,
    extract_features_from_variant,
    ACMG_GENES
)

__all__ = [
    "PathogenicityPredictor",
    "VariantFeatures", 
    "extract_features_from_variant",
    "ACMG_GENES"
]
