# data_filtering.py
"""
Модуль фильтрации данных для GenoScope.
Содержит универсальные фильтры по условиям, выбросам, корреляции, кластеризации и пайплайн-фильтрацию.
"""

import pandas as pd
import numpy as np
import traceback

from sklearn.cluster import KMeans
from utils.outlier_utils import get_outlier_mask

###################################
# I. Базовые методы
###################################


def filter_by_multiple_conditions(
    df: pd.DataFrame, conditions: list[str]
) -> pd.DataFrame:
    """
    Фильтрация строк на основе нескольких условий (df.query).

    Пример:
        conditions=["A > 10","B < 5"] → df.query("A > 10 & B < 5")

    Параметры:
        df (pd.DataFrame): Исходные данные.
        conditions (list of str): Список условий (строк для query).

    Возвращает:
        pd.DataFrame: Отфильтрованный DataFrame.
    """
    try:
        query_str = " & ".join(conditions)
        return df.query(query_str)
    except Exception as e:
        print(f"Ошибка фильтрации по условиям: {e}")
        traceback.print_exc()
        return df


def filter_by_custom_function(df: pd.DataFrame, func: callable) -> pd.DataFrame:
    """
    Фильтрация строк на основе пользовательской функции.

    Параметры:
        df (pd.DataFrame): Исходные данные.
        func (callable): Функция, принимающая строку (Series), возвращающая True/False.

    Возвращает:
        pd.DataFrame: Отфильтрованный DataFrame.
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


def filter_by_percentile(
    df: pd.DataFrame, column: str, lower: int = 5, upper: int = 95
) -> pd.DataFrame:
    """
    Удаляет строки, где значения столбца выходят за пределы заданных перцентилей.

    Параметры:
        df (pd.DataFrame): Исходные данные.
        column (str): Имя столбца.
        lower (float): Нижний перцентиль.
        upper (float): Верхний перцентиль.

    Возвращает:
        pd.DataFrame: Отфильтрованный DataFrame.
    """
    try:
        lb = np.percentile(df[column], lower)
        ub = np.percentile(df[column], upper)
        return df[(df[column] >= lb) & (df[column] <= ub)]
    except Exception as e:
        print(f"Error filtering by percentiles: {e}")
        traceback.print_exc()
        return df


def filter_outliers(
    df: pd.DataFrame,
    column: str | None = None,
    method: str = "iqr",
    threshold: float = 1.5,
) -> pd.DataFrame:
    """
    Удаляет выбросы из DataFrame на основе выбранного алгоритма.

    Параметры:
        df (pd.DataFrame): Данные для фильтрации.
        column (str или None): Столбец для фильтрации или None (по всем числовым столбцам).
        method (str): Метод обнаружения выбросов ('iqr', 'z-score', 'mahalanobis', 'isolation_forest').
        threshold (float): Параметр метода.

    Возвращает:
        pd.DataFrame: Очищенный DataFrame без выбросов.
    """
    mask = get_outlier_mask(df, column, method, threshold)
    if column:
        return df[~mask]
    else:
        if isinstance(mask, pd.DataFrame):
            return df[~mask.any(axis=1)]
        else:
            return df[~mask]


###################################
# III. Прочие универсальные методы
###################################


def filter_by_value_range(
    df: pd.DataFrame,
    column: str,
    min_value: float | None = None,
    max_value: float | None = None,
) -> pd.DataFrame:
    """
    Удаляет строки, где значение столбца выходит за пределы [min_value, max_value].

    Параметры:
        df (pd.DataFrame): Исходные данные.
        column (str): Имя столбца.
        min_value (float, optional): Нижняя граница (если None, не используется).
        max_value (float, optional): Верхняя граница (если None, не используется).

    Возвращает:
        pd.DataFrame: Отфильтрованный DataFrame.
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


def filter_by_nan(
    df: pd.DataFrame, how: str = "any", axis: int = 0, subset: list[str] | None = None
) -> pd.DataFrame:
    """
    Удаляет строки или столбцы с NaN (pandas dropna).

    Параметры:
        df (pd.DataFrame): Данные.
        how (str): 'any' — если есть хотя бы один NaN; 'all' — если все NaN.
        axis (int): 0 — строки, 1 — столбцы.
        subset (list, optional): Список столбцов для анализа NaN.

    Возвращает:
        pd.DataFrame: DataFrame без NaN.
    """
    try:
        return df.dropna(axis=axis, how=how, subset=subset)
    except Exception as e:
        print(f"Error in filter_by_nan: {e}")
        traceback.print_exc()
        return df


def filter_duplicates(
    df: pd.DataFrame, keep: str = "first", subset: list[str] | None = None
) -> pd.DataFrame:
    """
    Удаляет дубликаты строк (pandas drop_duplicates).

    Параметры:
        df (pd.DataFrame): Исходные данные.
        keep (str): 'first', 'last' или False.
        subset (list, optional): Список столбцов для сравнения.

    Возвращает:
        pd.DataFrame: DataFrame без дубликатов.
    """
    try:
        return df.drop_duplicates(subset=subset, keep=keep)
    except Exception as e:
        print(f"Error in filter_duplicates: {e}")
        traceback.print_exc()
        return df


def filter_by_correlation(
    df: pd.DataFrame, corr_threshold: float = 0.95
) -> pd.DataFrame:
    """
    Удаляет сильно коррелированные столбцы (feature selection).

    Параметры:
        df (pd.DataFrame): Данные.
        corr_threshold (float): Порог корреляции, выше которого удалять столбцы.

    Возвращает:
        pd.DataFrame: DataFrame без коррелированных столбцов.
    """
    try:
        corr_matrix = df.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        to_drop = [
            column for column in upper.columns if any(upper[column] > corr_threshold)
        ]
        return df.drop(columns=to_drop)
    except Exception as e:
        print(f"Error in filter_by_correlation: {e}")
        traceback.print_exc()
        return df


def filter_by_variance(df: pd.DataFrame, var_threshold: float = 0.0) -> pd.DataFrame:
    """
    Удаляет признаки (столбцы) с дисперсией <= var_threshold.

    Параметры:
        df (pd.DataFrame): Данные.
        var_threshold (float): Порог дисперсии (0 — только нулевая).

    Возвращает:
        pd.DataFrame: DataFrame без "пустых" признаков.
    """
    try:
        variances = df.var(numeric_only=True)
        low_var_cols = variances[variances <= var_threshold].index
        return df.drop(columns=low_var_cols)
    except Exception as e:
        print(f"Error in filter_by_variance: {e}")
        traceback.print_exc()
        return df


def filter_by_cluster(
    df: pd.DataFrame,
    method: str = "kmeans",
    n_clusters: int = 2,
    drop_outliers: bool = True,
) -> pd.DataFrame:
    """
    Удаляет точки, не вписывающиеся в кластеры (напр. методом KMeans).

    Параметры:
        df (pd.DataFrame): Исходные данные.
        method (str): Метод кластеризации ('kmeans' поддерживается).
        n_clusters (int): Число кластеров.
        drop_outliers (bool): Удалять ли выбросы вдали от центроида.

    Возвращает:
        pd.DataFrame: Отфильтрованный DataFrame.
    """
    try:
        if method != "kmeans":
            print("[filter_by_cluster] only 'kmeans' supported now.")
            return df

        num_df = df.select_dtypes(include=[np.number]).dropna(axis=0)
        if num_df.empty:
            return df

        km = KMeans(n_clusters=n_clusters, random_state=0)
        labels = km.fit_predict(num_df)
        centroids = km.cluster_centers_

        distances = []
        for i, row in enumerate(num_df.values):
            cluster_id = labels[i]
            center = centroids[cluster_id]
            dist = np.linalg.norm(row - center)
            distances.append(dist)
        distances = np.array(distances)

        if drop_outliers:
            mean_dist = distances.mean()
            std_dist = distances.std()
            cutoff = mean_dist + 2 * std_dist
            good_mask = distances <= cutoff
            good_idx = num_df.index[good_mask]
            return df.loc[good_idx]
        else:
            return df

    except Exception as e:
        print(f"Error in filter_by_cluster: {e}")
        traceback.print_exc()
        return df


###################################
# IV. Единый apply_filter_pipeline
###################################


def apply_filter_pipeline(df: pd.DataFrame, pipeline: list[dict]) -> pd.DataFrame:
    """
    Последовательно применяет несколько шагов фильтрации к DataFrame.

    Каждый шаг задаётся словарём, например:
        {"type": "outliers", "column": "X", "method": "iqr", ...}

    Параметры:
        df (pd.DataFrame): Данные.
        pipeline (list of dict): Описание шагов фильтрации.

    Возвращает:
        pd.DataFrame: DataFrame после всех фильтров.
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
                print(
                    f"[apply_filter_pipeline] Неизвестный тип фильтрации: {step_type}"
                )

        except Exception as e:
            print(f"[apply_filter_pipeline] Ошибка при шаге '{step_type}': {e}")
            traceback.print_exc()
            # Можно return df, если хотим прервать
            # return df

    return df
