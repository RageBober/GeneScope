# data_analysis/__init__.py
"""
Единая точка экспорта публичных API GenoScope.
"""

# ― Аналитика ―
from .analysis_core import (
    analyze_data,
    correlation_analysis,
    extract_pca,
    generate_statistics,
    select_features,
)
from .data_cleaning import (
    detect_outliers,
    handle_missing_values,
    remove_duplicates,
)

# ― Функции фильтрации ―
from .data_filtering import filter_outliers  # ← ИСПРАВЛЕНО: правильный модуль

__all__: list[str] = [
    # cleaning
    "remove_duplicates",
    "handle_missing_values",
    "detect_outliers",
    # filtering
    "filter_outliers",
    # analytics
    "select_features",
    "extract_pca",
    "correlation_analysis",
    "generate_statistics",
    "analyze_data",
]
