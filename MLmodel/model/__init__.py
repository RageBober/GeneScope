"""
MLmodel.model
-------------

Минимальный baseline-модель, чтобы код не падал при импорте.
Сейчас это RandomForestRegressor; при желании подмените позже.
"""

from pathlib import Path
import joblib
from sklearn.ensemble import RandomForestRegressor

_MODEL_PATH = Path(__file__).with_suffix(".joblib")


def get_model():
    """Загрузить модель, если она сохранена, или создать новую."""
    if _MODEL_PATH.exists():
        return joblib.load(_MODEL_PATH)
    return RandomForestRegressor(
        n_estimators=100,
        random_state=0,
        n_jobs=-1,
    )


__all__ = ["get_model"]
