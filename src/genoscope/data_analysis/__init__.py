# data_analysis/__init__.py
"""
Единая точка экспорта публичных API GenoScope.
"""

# ― Аналитика ―
from .analysis_core import analyze_data
from .analysis_core import correlation_analysis
from .analysis_core import extract_pca
from .analysis_core import generate_statistics
from .analysis_core import select_features
from .data_cleaning import detect_outliers
from .data_cleaning import handle_missing_values
from .data_cleaning import remove_duplicates

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
