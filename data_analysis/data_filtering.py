from __future__ import annotations
import traceback
from typing import Final


import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.cluster import KMeans
from sklearn.covariance import EmpiricalCovariance
from sklearn.ensemble import IsolationForest

###################################
# I. Базовые методы
###################################


def filter_by_multiple_conditions(df, conditions):
    """
    Фильтрация строк на основе нескольких условий (df.query).
    Пример: conditions=["A > 10","B < 5"] → df.query("A > 10 & B < 5").
    """
    try:
        query_str = " & ".join(conditions)
        return df.query(query_str)
    except Exception as e:
        print(f"Ошибка фильтрации по условиям: {e}")
        traceback.print_exc()
        return df


def filter_by_custom_function(df, func):
    """
    Фильтрация строк на основе пользовательской логики.
    func(row) -> True/False, где row - строка DataFrame.
    """
    try:
        return df[df.apply(func, axis=1)]
    except Exception as e:
        print(f"Ошибка пользовательской фильтрации: {e}")
        traceback.print_exc()
        return df


###################################
# II. Методы, основанные на статистике
###################################


def filter_by_percentile(df, column, lower=5, upper=95):
    """
    Удаление строк, где значения столбца выходят за [lower_percentile, upper_percentile].
    """
    try:
        lb = np.percentile(df[column], lower)
        ub = np.percentile(df[column], upper)
        return df[(df[column] >= lb) & (df[column] <= ub)]
    except Exception as e:
        print(f"Error filtering by percentiles: {e}")
        traceback.print_exc()
        return df


# -----------------------------------------------------------------------------
# Функция-обёртка
# -----------------------------------------------------------------------------
def filter_outliers(
    df: pd.DataFrame,
    column: str | None = None,
    *,
    method: str = "iqr",
    threshold: float = 1.5,
) -> pd.DataFrame:
    """
    Удаляет выбросы (строки) по одному из четырёх методов.

    Parameters
    ----------
    df : pd.DataFrame
    column : str | None
        Для "iqr" и "z-score" — числовой столбец, по которому ищем выбросы.
        Для "mahalanobis" и "isolation_forest" игнорируется (берутся все
        числовые фичи).
    method : {"iqr", "z-score", "mahalanobis", "isolation_forest"}
    threshold : float
        ▸ IQR — множитель (обычно 1.5)  
        ▸ Z-score — максимальный |z| (обычно 3)  
        ▸ Mahalanobis — порог расстояния (≈ χ²-квантиль)  
        ▸ IsolationForest — contamination = threshold (доля выбросов)

    Returns
    -------
    pd.DataFrame
        df без строк-выбросов
    """
    method: Final = method.lower()

    # --- sanity checks -------------------------------------------------------
    if method in {"iqr", "z-score"} and not column:
        raise ValueError(
            f"[filter_outliers] Для метода “{method}” нужно указать column."
        )
    if method not in {"iqr", "z-score", "mahalanobis", "isolation_forest"}:
        raise ValueError(f"[filter_outliers] Unsupported method: {method}")

    try:
        # ------------------------------------------------------------------ #
        if method == "iqr":
            q1, q3 = df[column].quantile([0.25, 0.75])
            iqr = q3 - q1
            lb, ub = q1 - threshold * iqr, q3 + threshold * iqr
            return df[(df[column] >= lb) & (df[column] <= ub)]

        # ------------------------------------------------------------------ #
        if method == "z-score":
            z = zscore(df[column].astype(float), nan_policy="omit")
            return df[np.abs(z) <= threshold]

        # ------------------------------------------------------------------ #
        if method == "mahalanobis":
            num_df = df.select_dtypes(include=[np.number]).dropna(axis=0)
            if num_df.empty or len(num_df) < 3:
                # Махаланобис плохо работает на 1-2 строках
                return df
            cov = EmpiricalCovariance().fit(num_df)
            distances = cov.mahalanobis(num_df)
            mask_keep = distances <= threshold
            return df.loc[num_df.index[mask_keep]]

        # ------------------------------------------------------------------ #
        # isolation_forest
        cont = float(threshold) if 0 < threshold < 1 else 0.05
        num_df = df.select_dtypes(include=[np.number]).dropna(axis=0)
        if num_df.empty:
            return df
        iso = IsolationForest(contamination=cont, random_state=0)
        preds = iso.fit_predict(num_df)  # +1 — норма, -1 — аномалия
        mask_keep = preds == 1
        return df.loc[num_df.index[mask_keep]]

    except Exception as exc:  # <- ловим, чтобы не уронить пайплайн
        print(f"[filter_outliers] Error: {exc}")
        traceback.print_exc()
        return df



###################################
# III. Прочие универсальные методы
###################################


def filter_by_value_range(df, column, min_value=None, max_value=None):
    """
    Удаление строк, где значение столбца column выходит за диапазон [min_value, max_value].
    Если min_value=None, не ограничиваем снизу;
    если max_value=None, не ограничиваем сверху.
    """
    try:
        res = df
        if min_value is not None:
            res = res[res[column] >= min_value]
        if max_value is not None:
            res = res[res[column] <= max_value]
        return res
    except Exception as e:
        print(f"Error in filter_by_value_range: {e}")
        traceback.print_exc()
        return df


def filter_by_nan(df, how="any", axis=0, subset=None):
    """
    Удаление строк или столбцов с NaN (pandas dropna).
    ПАРАМЕТРЫ:
      how='any'/'all'
      axis=0 => строки, axis=1 => столбцы
      subset: список столбцов
    """
    try:
        return df.dropna(axis=axis, how=how, subset=subset)
    except Exception as e:
        print(f"Error in filter_by_nan: {e}")
        traceback.print_exc()
        return df


def filter_duplicates(df, keep="first", subset=None):
    """
    Удаление дубликатов строк (pandas drop_duplicates).
    keep='first'/'last'/False,
    subset: список столбцов
    """
    try:
        return df.drop_duplicates(subset=subset, keep=keep)
    except Exception as e:
        print(f"Error in filter_duplicates: {e}")
        traceback.print_exc()
        return df


def filter_by_correlation(df, corr_threshold=0.95):
    """
    Удаление сильно коррелированных столбцов (feature selection).
    Оставляем по одному из группы столбцов, корелл. выше corr_threshold.
    """
    try:
        corr_matrix = df.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        to_drop = [column for column in upper.columns if any(upper[column] > corr_threshold)]
        return df.drop(columns=to_drop)
    except Exception as e:
        print(f"Error in filter_by_correlation: {e}")
        traceback.print_exc()
        return df


def filter_by_variance(df, var_threshold=0.0):
    """
    Удаление признаков (столбцов) с дисперсией <= var_threshold.
    var_threshold=0 -> удаляем столбцы, где var=0.
    """
    try:
        variances = df.var(numeric_only=True)
        low_var_cols = variances[variances <= var_threshold].index
        return df.drop(columns=low_var_cols)
    except Exception as e:
        print(f"Error in filter_by_variance: {e}")
        traceback.print_exc()
        return df


def filter_by_cluster(df, method="kmeans", n_clusters=2, drop_outliers=True):
    """
    Удаление точек, которые "не вписываются" в кластеры.
    Здесь реальная реализация KMeans.
    ПАРАМЕТРЫ:
      method='kmeans'
      n_clusters=2
      drop_outliers=True (предположим, убираем точки, у которых расстояние > 2*std от центроида)
    """
    try:
        if method != "kmeans":
            print("[filter_by_cluster] only 'kmeans' supported now.")
            return df

        # Берём все числовые столбцы, убираем NaN
        num_df = df.select_dtypes(include=[np.number]).dropna(axis=0)
        if num_df.empty:
            return df

        km = KMeans(n_clusters=n_clusters, random_state=0)
        labels = km.fit_predict(num_df)
        centroids = km.cluster_centers_

        # Для каждой точки считаем расстояние до её центроида
        distances = []
        for i, row in enumerate(num_df.values):
            cluster_id = labels[i]
            center = centroids[cluster_id]
            dist = np.linalg.norm(row - center)
            distances.append(dist)
        distances = np.array(distances)

        # Примем, что "слишком далеко" = dist > mean+2*std
        # (очень грубо).
        if drop_outliers:
            mean_dist = distances.mean()
            std_dist = distances.std()
            cutoff = mean_dist + 2 * std_dist
            good_mask = distances <= cutoff
            good_idx = num_df.index[good_mask]
            return df.loc[good_idx]
        else:
            return df  # не убираем "дальних" — просто не фильтруем

    except Exception as e:
        print(f"Error in filter_by_cluster: {e}")
        traceback.print_exc()
        return df


###################################
# IV. Единый apply_filter_pipeline
###################################


def apply_filter_pipeline(df, pipeline):
    """
    Применяет несколько шагов фильтрации подряд к DataFrame.
    Каждый шаг — словарь вида:
      { "type": "...", ... }
    """
    for step in pipeline:
        step_type = step.get("type")
        try:
            if step_type == "multiple_conditions":
                conditions = step.get("conditions", [])
                df = filter_by_multiple_conditions(df, conditions)

            elif step_type == "custom_function":
                func = step.get("func")
                if func is not None:
                    df = filter_by_custom_function(df, func)

            elif step_type == "percentile":
                column = step.get("column")
                lower = step.get("lower", 5)
                upper = step.get("upper", 95)
                if column:
                    df = filter_by_percentile(df, column, lower, upper)

            elif step_type == "outliers":
                column = step.get("column")
                method = step.get("method", "iqr")
                threshold = step.get("threshold", 1.5)
                df = filter_outliers(df, column, method=method, threshold=threshold)

            elif step_type == "value_range":
                col = step.get("column")
                min_val = step.get("min_value", None)
                max_val = step.get("max_value", None)
                if col:
                    df = filter_by_value_range(df, col, min_val, max_val)

            elif step_type == "nan_filter":
                how = step.get("how", "any")
                axis = step.get("axis", 0)
                subset = step.get("subset", None)
                df = filter_by_nan(df, how=how, axis=axis, subset=subset)

            elif step_type == "duplicates_filter":
                keep = step.get("keep", "first")
                subset = step.get("subset", None)
                df = filter_duplicates(df, keep=keep, subset=subset)

            elif step_type == "correlation_filter":
                thr = step.get("corr_threshold", 0.95)
                df = filter_by_correlation(df, corr_threshold=thr)

            elif step_type == "variance_filter":
                var_thr = step.get("var_threshold", 0.0)
                df = filter_by_variance(df, var_threshold=var_thr)

            elif step_type == "cluster_filter":
                method = step.get("method", "kmeans")
                n_clust = step.get("n_clusters", 2)
                drop_out = step.get("drop_outliers", True)
                df = filter_by_cluster(
                    df, method=method, n_clusters=n_clust, drop_outliers=drop_out
                )

            else:
                print(f"[apply_filter_pipeline] Неизвестный тип фильтрации: {step_type}")

        except Exception as e:
            print(f"[apply_filter_pipeline] Ошибка при шаге '{step_type}': {e}")
            traceback.print_exc()
            # Можно return df, если хотим прервать
            # return df

    return df
