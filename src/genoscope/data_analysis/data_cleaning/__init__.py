# src/data_analysis/data_cleaning/__init__.py

from __future__ import annotations

from .data_cleaning_ml import detect_outliers as _detect_outliers_ml  # ML-блок
from .data_cleaning_ml import fill_missing_with_ml  # ML-блок
from .data_cleaning_statistic import detect_outliers as _detect_outliers_stat  # статистика
from .data_cleaning_statistic import handle_missing_values  # статистика
from .data_cleaning_statistic import remove_duplicates  # статистика


# ----- единая точка входа -------------------------------------------------
def detect_outliers(
    df,
    column: str | None = None,
    *,
    method: str = "iqr",  # type: ignore[arg-type]
    threshold: float = 1.5,
):
    """
    Обёртка, автоматически перенаправляющая вызов:
    • IQR / Z-score / Mahalanobis  → статистический модуль;
    • Isolation-Forest             → ML-модуль.
    """
    if method.lower() == "isolation_forest":
        return _detect_outliers_ml(
            df, column, method="isolation_forest", threshold=threshold
        )
    return _detect_outliers_stat(df, column, method=method, threshold=threshold)  # type: ignore[arg-type]


__all__ = [
    # статистика
    "remove_duplicates",
    "handle_missing_values",
    # ML
    "fill_missing_with_ml",
    # универсальная функция
    "detect_outliers",
]
