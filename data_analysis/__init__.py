# data_analysis/__init__.py

# --- Загрузка / приём данных -------------------------------------------------
from .data_ingestion import load_data
from .data_filtering import filter_outliers

from .data_cleaning import (
    remove_duplicates,
    handle_missing_values,
    detect_outliers,
)

# --- Фильтрация выбросов -----------------------------------------------------
#  filter_outliers живёт в data_filtering.py, импортируем оттуда
# -------------------------------------------------------------------
from .analysis_core import extract_pca, select_features, correlation_analysis
from .visualization import plot_pca
from .sequence_analysis import encode_sequences

__all__ = [
    "load_data",
    "remove_duplicates",
    "handle_missing_values",
    "detect_outliers",
    "filter_outliers",
    "extract_pca",
    "select_features",
    "correlation_analysis",
    "plot_pca",
    "encode_sequences",
]
