"""
data_filters.pipeline
---------------------

Единая точка входа для последовательного применения низкоуровневых
фильтров к DataFrame.

Пример:
    from data_filters.pipeline import run_filters

    df_clean = run_filters(
        df,
        steps=[
            {"name": "expression", "params": {"min_expr": 5}},
            {"name": "variants",   "params": {"impact": "HIGH"}},
            {"name": "outliers",   "params": {"column": "depth", "method": "iqr"}},
        ]
    )
"""

from typing import List, Dict, Any
import pandas as pd

# импортируем конкретные фильтры
from .filter_expression import filter_expression
from .filter_meta_epi import filter_meta_epi
from .filter_ngs import filter_ngs
from .filter_variants import filter_variants
from ..data_filtering import filter_outliers

# Регистрируем доступные фильтры в словаре
_FILTER_REGISTRY = {
    "expression": filter_expression,
    "meta_epi": filter_meta_epi,
    "ngs": filter_ngs,
    "variants": filter_variants,
    "outliers": filter_outliers,
}


def run_filters(df: pd.DataFrame, steps: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Последовательно применяет фильтры.

    Parameters
    ----------
    df : pd.DataFrame
        Входные данные.
    steps : list[dict]
        Каждый элемент: {"name": <ключ из _FILTER_REGISTRY>, "params": {...}}

    Returns
    -------
    pd.DataFrame
        Данные после всех фильтров.
    """
    result = df.copy()

    for step in steps:
        name = step["name"]
        params = step.get("params", {})

        fn = _FILTER_REGISTRY.get(name)
        if fn is None:
            raise ValueError(f"[pipeline] Неизвестный фильтр: '{name}'")

        result = fn(result, **params)

    return result


__all__ = ["run_filters"]
