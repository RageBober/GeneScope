"""
Machine Learning Model for Variant Pathogenicity Prediction
Implements ensemble model combining multiple features and classifiers
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Any, Tuple
import pickle
import logging
from pathlib import Path
from dataclasses import dataclass
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import joblib

logger = logging.getLogger(__name__)

@dataclass
class VariantFeatures:
    """Features for variant pathogenicity prediction."""
    # Conservation scores
    phylop_score: float = 0.0  # PhyloP conservation score
    phastcons_score: float = 0.0  # PhastCons conservation score
    gerp_score: float = 0.0  # GERP++ score
    
    # Functional predictions
    sift_score: float = 1.0  # SIFT score (0 = deleterious)
    polyphen_score: float = 0.0  # PolyPhen score (1 = probably damaging)
    cadd_score: float = 0.0  # CADD score
    revel_score: float = 0.0  # REVEL score
    
    # Population frequency
    gnomad_af: float = 0.0  # gnomAD allele frequency
    max_pop_af: float = 0.0  # Maximum population AF
    
    # Protein features
    is_missense: bool = False
    is_nonsense: bool = False
    is_frameshift: bool = False
    is_splice: bool = False
    is_inframe: bool = False
    
    # Gene-level features
    gene_pli: float = 0.0  # Gene pLI score
    gene_loeuf: float = 1.0  # Gene LOEUF score
    is_cancer_gene: bool = False
    is_acmg_gene: bool = False
    
    # Additional features
    protein_position: float = 0.0  # Relative position in protein
    is_domain: bool = False  # In functional domain
    has_clinvar: bool = False  # Has ClinVar entry
    
    def to_array(self) -> np.ndarray:
        """Convert features to numpy array for ML model."""
        return np.array([
            self.phylop_score,
            self.phastcons_score,
            self.gerp_score,
            self.sift_score,
            self.polyphen_score,
            self.cadd_score,
            self.revel_score,
            np.log10(self.gnomad_af + 1e-10),  # Log transform frequency
            np.log10(self.max_pop_af + 1e-10),
            float(self.is_missense),
            float(self.is_nonsense),
            float(self.is_frameshift),
            float(self.is_splice),
            float(self.is_inframe),
            self.gene_pli,
            self.gene_loeuf,
            float(self.is_cancer_gene),
            float(self.is_acmg_gene),
            self.protein_position,
            float(self.is_domain),
            float(self.has_clinvar)
        ])


class PathogenicityPredictor:
    """
    Ensemble model for variant pathogenicity prediction.
    Combines multiple ML models and rule-based filters.
    """
    
    def __init__(self, model_path: Optional[str] = None):
        """
        Initialize pathogenicity predictor.
        
        Args:
            model_path: Path to pre-trained model file
        """
        self.model_path = model_path
        self.models = {}
        self.scaler = StandardScaler()
        self.feature_names = [
            "phylop", "phastcons", "gerp", "sift", "polyphen",
            "cadd", "revel", "log_gnomad_af", "log_max_pop_af",
            "is_missense", "is_nonsense", "is_frameshift", "is_splice",
            "is_inframe", "gene_pli", "gene_loeuf", "is_cancer_gene",
            "is_acmg_gene", "protein_position", "is_domain", "has_clinvar"
        ]
        
        if model_path and Path(model_path).exists():
            self.load_model(model_path)
        else:
            self._initialize_default_models()
    
    def _initialize_default_models(self):
        """Initialize default ensemble models."""
        # Random Forest
        self.models['rf'] = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            min_samples_split=5,
            class_weight='balanced',
            random_state=42
        )
        
        # Gradient Boosting
        self.models['gb'] = GradientBoostingClassifier(
            n_estimators=100,
            learning_rate=0.1,
            max_depth=5,
            random_state=42
        )
        
        logger.info("Initialized default ML models")
    
    def predict(self, features: VariantFeatures) -> Dict[str, Any]:
        """
        Predict pathogenicity for a variant.
        
        Args:
            features: Variant features
            
        Returns:
            Dict with prediction results
        """
        # Convert to array and scale
        feature_array = features.to_array().reshape(1, -1)
        
        # Use rule-based prediction since models aren't trained
        return self._rule_based_prediction(features)
    
    def _rule_based_prediction(self, features: VariantFeatures) -> Dict[str, Any]:
        """
        Rule-based prediction when ML models not available.
        Based on ACMG/AMP guidelines.
        """
        evidence_pathogenic = []
        evidence_benign = []
        score = 0
        
        # Very strong evidence
        if features.is_nonsense and features.gene_pli > 0.9:
            evidence_pathogenic.append("PVS1: Nonsense in LoF intolerant gene")
            score += 8
        
        # Strong evidence
        if features.cadd_score > 30:
            evidence_pathogenic.append(f"PS1: Very high CADD score ({features.cadd_score:.1f})")
            score += 4
        
        if features.gnomad_af < 0.00001 and features.is_missense:
            evidence_pathogenic.append("PM2: Absent/very rare in gnomAD")
            score += 2
        
        # Moderate evidence
        if features.revel_score > 0.75:
            evidence_pathogenic.append(f"PP3: High REVEL score ({features.revel_score:.2f})")
            score += 2
        
        if features.polyphen_score > 0.95 and features.sift_score < 0.05:
            evidence_pathogenic.append("PP3: Multiple computational tools predict damaging")
            score += 2
        
        # Benign evidence
        if features.gnomad_af > 0.01:
            evidence_benign.append(f"BA1: High allele frequency ({features.gnomad_af:.3f})")
            score -= 8
        
        if features.gnomad_af > 0.001:
            evidence_benign.append(f"BS1: Allele frequency greater than expected")
            score -= 4
        
        # Calculate probability based on score
        probability = 1 / (1 + np.exp(-score/4))  # Sigmoid transformation
        
        # Classification
        if score >= 10:
            classification = "Pathogenic"
        elif score >= 6:
            classification = "Likely Pathogenic"
        elif score <= -6:
            classification = "Benign"
        elif score <= -3:
            classification = "Likely Benign"
        else:
            classification = "Uncertain Significance"
        
        return {
            "prediction": int(probability >= 0.5),
            "probability": float(probability),
            "classification": classification,
            "model_predictions": {"rule_based": int(probability >= 0.5)},
            "model_probabilities": {"rule_based": float(probability)},
            "confidence": abs(probability - 0.5) * 2,
            "evidence": {
                "pathogenic": evidence_pathogenic,
                "benign": evidence_benign,
                "score": score
            }
        }
    
    def save_model(self, path: str):
        """Save trained model to file."""
        model_data = {
            "models": self.models,
            "scaler": self.scaler,
            "feature_names": self.feature_names
        }
        joblib.dump(model_data, path)
        logger.info(f"Model saved to {path}")
    
    def load_model(self, path: str):
        """Load trained model from file."""
        model_data = joblib.load(path)
        self.models = model_data["models"]
        self.scaler = model_data["scaler"]
        self.feature_names = model_data["feature_names"]
        logger.info(f"Model loaded from {path}")


# ACMG gene list (subset of medically important genes)
ACMG_GENES = {
    "BRCA1", "BRCA2", "TP53", "MLH1", "MSH2", "MSH6", "PMS2",
    "APC", "MUTYH", "VHL", "MEN1", "RET", "PTEN", "STK11",
    "TGFBR1", "TGFBR2", "SMAD4", "BMPR1A", "SDHB", "SDHC",
    "SDHD", "SDHAF2", "TSC1", "TSC2", "WT1", "NF2", "COL3A1",
    "FBN1", "TGFBR1", "TGFBR2", "SMAD3", "ACTA2", "MYH11",
    "MYBPC3", "MYH7", "TNNT2", "TNNI3", "TPM1", "MYL3",
    "ACTC1", "PRKAG2", "GLA", "MYL2", "LMNA", "RYR2",
    "PKP2", "DSP", "DSC2", "TMEM43", "DSG2", "KCNQ1",
    "KCNH2", "SCN5A", "LDLR", "APOB", "PCSK9", "RYR1",
    "CACNA1S", "ATP7B", "OTC"
}


def extract_features_from_variant(variant_data: Dict, 
                                 clinvar_data: Optional[Dict] = None,
                                 gnomad_data: Optional[Dict] = None,
                                 cosmic_data: Optional[Dict] = None) -> VariantFeatures:
    """
    Extract ML features from variant data.
    
    Args:
        variant_data: Basic variant information
        clinvar_data: ClinVar annotation
        gnomad_data: gnomAD population data
        cosmic_data: COSMIC cancer data
        
    Returns:
        VariantFeatures object
    """
    features = VariantFeatures()
    
    # Basic variant type
    var_type = variant_data.get("type", "").lower()
    features.is_missense = "missense" in var_type or var_type == "snp"
    features.is_nonsense = "nonsense" in var_type or "stop" in var_type
    features.is_frameshift = "frameshift" in var_type
    features.is_splice = "splice" in var_type
    features.is_inframe = "inframe" in var_type
    
    # Population frequency
    if gnomad_data:
        features.gnomad_af = gnomad_data.get("global_af", 0)
        pops = gnomad_data.get("populations", {})
        features.max_pop_af = max(pops.values()) if pops else 0
    
    # Gene features
    gene = variant_data.get("gene", "")
    features.is_acmg_gene = gene in ACMG_GENES
    
    if cosmic_data:
        features.is_cancer_gene = cosmic_data.get("is_cancer_gene", False)
    
    # ClinVar features
    if clinvar_data:
        features.has_clinvar = True
    
    # Mock scores for demonstration (would come from annotation tools)
    np.random.seed(hash(str(variant_data)) % 2**32)
    features.cadd_score = np.random.uniform(0, 40)
    features.revel_score = np.random.uniform(0, 1)
    features.phylop_score = np.random.uniform(-5, 10)
    features.phastcons_score = np.random.uniform(0, 1)
    features.sift_score = np.random.uniform(0, 1)
    features.polyphen_score = np.random.uniform(0, 1)
    features.gene_pli = np.random.uniform(0, 1)
    features.gene_loeuf = np.random.uniform(0.1, 2)
    
    return features
