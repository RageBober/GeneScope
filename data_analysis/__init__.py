# data_analysis/__init__.py
"""
Единая точка экспорта публичных API GenoScope.
"""

from .data_cleaning import (
    remove_duplicates,
    handle_missing_values,
    detect_outliers,
)

# ― Функции фильтрации ―
from .data_filtering import filter_outliers  # ← ИСПРАВЛЕНО: правильный модуль

# ― Аналитика ―
from .analysis_core import (
    select_features,
    extract_pca,
    correlation_analysis,
    generate_statistics,
    analyze_data,
)

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
