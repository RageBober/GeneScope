# filter_expression.py

import numpy as np
import pandas as pd

def filter_by_expression_level(df, column="TPM", threshold=1.0, how="above"):
    """
    Удаляем гены (строки), у которых средний уровень экспрессии (в колонке column) 
    ниже (или выше) заданного порога threshold.

    Параметры:
    - df (pd.DataFrame): DataFrame, где строки — гены, столбцы — образцы (либо наоборот)
    - column (str): название столбца, где хранится TPM / CPM
    - threshold (float): пороговое значение
    - how (str): 'above' (оставляем только >= threshold) или 'below' (оставляем только <= threshold)

    Возвращает:
    - pd.DataFrame (отфильтрованный)
    """
    if column not in df.columns:
        print(f"[filter_by_expression_level] Колонка {column} не найдена, пропускаем фильтр.")
        return df

    if how == "above":
        return df[df[column] >= threshold]
    elif how == "below":
        return df[df[column] <= threshold]
    else:
        print(f"[filter_by_expression_level] Неизвестный параметр how: {how}")
        return df


def filter_by_zero_count(df, columns=None, max_zero_fraction=0.5):
    """
    Удаляем гены, у которых доля нулевых значений > max_zero_fraction
    Например, если max_zero_fraction=0.5, то если 50% образцов = 0, ген удаляется.

    Параметры:
    - df (pd.DataFrame): строки — гены, столбцы — образцы
    - columns (list|None): список столбцов, среди которых искать нули. 
                           Если None — берём все числовые столбцы.
    - max_zero_fraction (float): допустимая доля нулей

    Возвращает:
    - pd.DataFrame
    """
    if columns is None:
        # Берём все столбцы float/int
        columns = df.select_dtypes(include=[np.number]).columns

    mask = []
    for idx, row in df.iterrows():
        total = len(columns)
        zero_count = (row[columns] == 0).sum()
        fraction = zero_count / total if total > 0 else 0
        # Если fraction <= max_zero_fraction, то оставляем
        if fraction <= max_zero_fraction:
            mask.append(True)
        else:
            mask.append(False)

    return df[mask]


def filter_by_fold_change(df, column_fc="log2FC", min_fc=1.0):
    """
    Отбрасываем строки, у которых |log2FC| < min_fc.
    То есть если log2FC=1 -> это кратное изменение в 2 раза.

    Параметры:
    - df (pd.DataFrame): 
    - column_fc (str): столбец с log2FC
    - min_fc (float): минимальный порог
    Возвращает:
    - pd.DataFrame
    """
    if column_fc not in df.columns:
        print(f"[filter_by_fold_change] Колонка {column_fc} не найдена, пропускаем.")
        return df

    return df[np.abs(df[column_fc]) >= min_fc]


def filter_by_significance(df, column_p="pvalue", alpha=0.05, method="pvalue"):
    """
    Фильтрация по статистической значимости:
    - method="pvalue": оставляем строки, где pvalue <= alpha
    - method="FDR"   : возможно, есть столбец 'padj' или 'FDR'

    Параметры:
    - df
    - column_p: название столбца с p-value или FDR
    - alpha: порог
    - method: 'pvalue' или 'FDR'
    """
    if column_p not in df.columns:
        print(f"[filter_by_significance] Колонка {column_p} не найдена, пропускаем.")
        return df

    if method.lower() == "pvalue":
        return df[df[column_p] <= alpha]
    elif method.lower() == "fdr":
        return df[df[column_p] <= alpha]
    else:
        print(f"[filter_by_significance] Неизвестный method: {method}.")
        return df


def filter_by_biomarker_bases(df, gene_column="gene_id", biomarker_list=None):
    """
    Оставить в DataFrame только те гены, которые есть в списке биомаркеров/путей.
    Например, GO, KEGG, Reactome.
    Параметры:
    - df
    - gene_column: как называется столбец с идентификатором гена
    - biomarker_list: список (set) генов
    """
    if biomarker_list is None:
        # Нечего фильтровать
        return df

    if gene_column not in df.columns:
        print(f"[filter_by_biomarker_bases] Колонка {gene_column} не найдена, пропускаем.")
        return df

    return df[df[gene_column].isin(biomarker_list)]

'''о фильтрации 
Стоит ли сейчас допиливать этот модуль дальше
Если вы уже знаете, что переходите к pipeline/UI, в котором пользователь сам выберет «метод=AVG, threshold=1.0, require≥2 samples» – да, имеет смысл.

Если всё ещё в процессе «глобальной» оптимизации (Dask/C++), можно сначала набросать методы «как есть», а потом, при необходимости, переписать на Dask.'''