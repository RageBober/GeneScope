from __future__ import annotations

import logging
from collections.abc import Callable
from collections.abc import Sequence
from typing import Any
from typing import Literal
from typing import TypeAlias

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.stats import zscore
from sklearn.cluster import KMeans
from sklearn.covariance import EmpiricalCovariance
from sklearn.ensemble import IsolationForest

from genoscope.data_analysis._typing import BoolS

logger = logging.getLogger(__name__)

DF: TypeAlias = pd.DataFrame


###############################################################################
# I. Базовые методы
###############################################################################


def filter_by_multiple_conditions(df: DF, conditions: Sequence[str]) -> DF:
    """
    Фильтрация строк по списку условий (pandas.query).
    """
    try:
        query_str = " & ".join(conditions)
        return df.query(query_str)
    except Exception as exc:
        logger.error("filter_by_multiple_conditions: %s", exc, exc_info=True)
        return df


def filter_by_custom_function(df: DF, func: Callable[[pd.Series[Any]], bool]) -> DF:
    """
    Фильтрация строк по пользовательской функции: func(row) -> bool.
    """
    try:
        mask: BoolS = df.apply(func, axis=1)
        return df.loc[mask]
    except Exception as exc:
        logger.error("filter_by_custom_function: %s", exc, exc_info=True)
        return df


###############################################################################
# II. Методы, основанные на статистике
###############################################################################


def filter_by_percentile(
    df: DF,
    column: str,
    *,
    lower: float = 5.0,
    upper: float = 95.0,
) -> DF:
    """
    Удаление строк с значением столбца вне [lower, upper] percentiles.
    """
    try:
        values = pd.to_numeric(df[column], errors="coerce")
        percentiles = np.nanpercentile(values, [lower, upper])
        if np.isscalar(percentiles):
            # Handle edge case where all values are NaN
            return df
        lb, ub = percentiles[0], percentiles[1]
        mask: BoolS = (values >= lb) & (values <= ub)
        return df.loc[mask]
    except Exception as exc:
        logger.error("filter_by_percentile: %s", exc, exc_info=True)
        return df


Method = Literal["iqr", "z-score", "mahalanobis", "isolation_forest"]


def filter_outliers(
    df: DF,
    column: str | None = None,
    *,
    method: Method = "iqr",
    threshold: float = 1.5,
    contamination: float = 0.05,
) -> DF:
    """
    Удаление выбросов по одному из методов: IQR, Z-score, Mahalanobis, IsolationForest.
    Для IQR/Z-score требуется column.
    """
    algo = method.lower()
    if algo in {"iqr", "z-score"} and column is None:
        raise ValueError(f"filter_outliers: column required for {algo}")

    try:
        if algo == "iqr":
            arr = pd.to_numeric(df[column], errors="coerce").dropna()
            q1, q3 = arr.quantile([0.25, 0.75])
            iqr_val = float(q3 - q1)
            lb, ub = q1 - threshold * iqr_val, q3 + threshold * iqr_val
            mask_iqr: BoolS = (df[column] >= lb) & (df[column] <= ub)
            return df.loc[mask_iqr]

        if algo == "z-score":
            arr = pd.to_numeric(df[column], errors="coerce").astype(float)
            z_arr: NDArray[np.floating[Any]] = zscore(arr, nan_policy="omit")
            mask_z: BoolS = pd.Series(np.abs(z_arr) <= threshold, index=df.index)
            return df.loc[mask_z]

        if algo == "mahalanobis":
            num_df = df.select_dtypes(include=[np.number]).dropna()
            if num_df.shape[0] < 3:
                return df
            dist = EmpiricalCovariance().fit(num_df).mahalanobis(num_df)
            mask_mh: BoolS = pd.Series(dist <= threshold, index=num_df.index)
            return df.loc[mask_mh]

        # IsolationForest
        num_df = df.select_dtypes(include=[np.number]).dropna()
        if num_df.empty:
            return df
        iso = IsolationForest(contamination=contamination, random_state=0)
        preds = iso.fit_predict(num_df)
        mask_if: BoolS = pd.Series(preds == 1, index=num_df.index)
        return df.loc[mask_if]

    except Exception as exc:
        logger.error("filter_outliers: %s", exc, exc_info=True)
        return df


###############################################################################
# III. Универсальные методы
###############################################################################


def filter_by_value_range(
    df: DF,
    column: str,
    min_value: float | None = None,
    max_value: float | None = None,
) -> DF:
    """
    Оставляет строки, где column в диапазоне [min_value, max_value].
    """
    try:
        res = df
        if min_value is not None:
            res = res.loc[res[column] >= min_value]
        if max_value is not None:
            res = res.loc[res[column] <= max_value]
        return res
    except Exception as exc:
        logger.error("filter_by_value_range: %s", exc, exc_info=True)
        return df


def filter_by_nan(
    df: DF,
    *,
    how: Literal["any", "all"] = "any",
    axis: Literal[0, 1] = 0,
    subset: Sequence[str] | None = None,
) -> DF:
    """
    Удаление NaN по параметрам pandas.dropna.
    """
    try:
        return df.dropna(axis=axis, how=how, subset=subset)
    except Exception as exc:
        logger.error("filter_by_nan: %s", exc, exc_info=True)
        return df


KeepLit = Literal["first", "last", False]


def filter_duplicates(
    df: DF,
    *,
    keep: KeepLit = "first",
    subset: Sequence[str] | None = None,
) -> DF:
    """
    Удаление дубликатов pandas.drop_duplicates.
    """
    try:
        return df.drop_duplicates(subset=subset, keep=keep)
    except Exception as exc:
        logger.error("filter_duplicates: %s", exc, exc_info=True)
        return df


def filter_by_correlation(df: DF, *, corr_threshold: float = 0.95) -> DF:
    """
    Удаление сильно скоррелированных столбцов.
    """
    try:
        corr = df.corr(numeric_only=True).abs()
        upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
        drop_cols = [c for c in upper.columns if any(upper[c] > corr_threshold)]
        return df.drop(columns=drop_cols)
    except Exception as exc:
        logger.error("filter_by_correlation: %s", exc, exc_info=True)
        return df


def filter_by_variance(df: DF, *, var_threshold: float = 0.0) -> DF:
    """
    Удаление столбцов с дисперсией <= var_threshold.
    """
    try:
        low_var = df.var(numeric_only=True) <= var_threshold
        return df.drop(columns=low_var[low_var].index)
    except Exception as exc:
        logger.error("filter_by_variance: %s", exc, exc_info=True)
        return df


def filter_by_cluster(
    df: DF,
    *,
    method: Literal["kmeans"] = "kmeans",
    n_clusters: int = 2,
    drop_outliers: bool = True,
) -> DF:
    """
    Кластеризация KMeans + удаление точек, далеких от центроидов.
    """
    if method.lower() != "kmeans":
        logger.warning("filter_by_cluster: только 'kmeans' поддерживается")
        return df
    try:
        num_df = df.select_dtypes(include=[np.number]).dropna()
        if num_df.empty:
            return df
        km = KMeans(n_clusters=n_clusters, random_state=0).fit(num_df)
        dist = np.linalg.norm(num_df.values - km.cluster_centers_[km.labels_], axis=1)
        if not drop_outliers:
            return df
        mask: BoolS = pd.Series(
            dist <= dist.mean() + 2 * dist.std(), index=num_df.index
        )
        return df.loc[mask]
    except Exception as exc:
        logger.error("filter_by_cluster: %s", exc, exc_info=True)
        return df


###############################################################################
# IV. Пайплайн
###############################################################################

PipelineStep = dict[str, Any]
Pipeline = Sequence[PipelineStep]


def apply_filter_pipeline(df: DF, pipeline: Pipeline) -> DF:
    """
    Применение фильтров по списку шагов: {'type': ..., ...}.
    """
    for step in pipeline:
        stype = step.get("type")
        if stype not in {
            "multiple_conditions",
            "custom_function",
            "percentile",
            "outliers",
            "value_range",
            "nan_filter",
            "duplicates_filter",
            "correlation_filter",
            "variance_filter",
            "cluster_filter",
        }:
            logger.warning("apply_filter_pipeline: unknown type '%s'", stype)
            continue
        try:
            if stype == "multiple_conditions":
                df = filter_by_multiple_conditions(df, step.get("conditions", []))
            elif stype == "custom_function":
                func = step.get("func")
                if func:
                    df = filter_by_custom_function(df, func)
            elif stype == "percentile":
                df = filter_by_percentile(
                    df,
                    column=step["column"],
                    lower=step.get("lower", 5.0),
                    upper=step.get("upper", 95.0),
                )
            elif stype == "outliers":
                df = filter_outliers(
                    df,
                    column=step.get("column"),
                    method=step.get("method", "iqr"),
                    threshold=step.get("threshold", 1.5),
                    contamination=step.get("contamination", 0.05),
                )
            elif stype == "value_range":
                col = step.get("column")
                if not isinstance(col, str):
                    logger.warning(
                        "apply_filter_pipeline: 'value_range' missing 'column'"
                    )
                    continue
                df = filter_by_value_range(
                    df,
                    column=col,
                    min_value=step.get("min_value"),
                    max_value=step.get("max_value"),
                )
            elif stype == "nan_filter":
                df = filter_by_nan(
                    df,
                    how=step.get("how", "any"),
                    axis=step.get("axis", 0),
                    subset=step.get("subset"),
                )
            elif stype == "duplicates_filter":
                df = filter_duplicates(
                    df,
                    keep=step.get("keep", "first"),
                    subset=step.get("subset"),
                )
            elif stype == "correlation_filter":
                df = filter_by_correlation(
                    df, corr_threshold=step.get("corr_threshold", 0.95)
                )
            elif stype == "variance_filter":
                df = filter_by_variance(
                    df, var_threshold=step.get("var_threshold", 0.0)
                )
            elif stype == "cluster_filter":
                df = filter_by_cluster(
                    df,
                    method=step.get("method", "kmeans"),
                    n_clusters=step.get("n_clusters", 2),
                    drop_outliers=step.get("drop_outliers", True),
                )
        except Exception as exc:
            logger.error("apply_filter_pipeline (%s): %s", stype, exc, exc_info=True)
    return df
