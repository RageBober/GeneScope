# src/data_analysis/data_cleaning/__init__.py

from __future__ import annotations

from .data_cleaning_statistic import (       # статистика
    remove_duplicates,
    handle_missing_values,
    detect_outliers as _detect_outliers_stat,
)
from .data_cleaning_ml import (              # ML-блок
    fill_missing_with_ml,
    detect_outliers as _detect_outliers_ml,
)

# ----- единая точка входа -------------------------------------------------
def detect_outliers(
    df,
    column: str | None = None,
    *,
    method: str = "iqr",
    threshold: float = 1.5,
):
    """
    Обёртка, автоматически перенаправляющая вызов:
    • IQR / Z-score / Mahalanobis  → статистический модуль;
    • Isolation-Forest             → ML-модуль.
    """
    if method.lower() == "isolation_forest":
        return _detect_outliers_ml(df, column, method="isolation_forest", threshold=threshold)
    return _detect_outliers_stat(df, column, method=method, threshold=threshold)

__all__ = [
    # статистика
    "remove_duplicates",
    "handle_missing_values",
    # ML
    "fill_missing_with_ml",
    # универсальная функция
    "detect_outliers",
]
