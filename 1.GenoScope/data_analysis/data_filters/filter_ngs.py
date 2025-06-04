# filter_ngs.py

import pandas as pd
import numpy as np

def filter_by_phred_quality(df, column="mean_phred", min_phred=20):
    """
    Удаляем риды, где среднее Phred-качество < min_phred.
    
    Параметры:
      df (DataFrame): содержащий столбец 'column'.
      column (str): название столбца, где хранится оценка качества (число).
      min_phred (int): порог, ниже которого риды считаем некачественными.
    
    Возвращает:
      DataFrame: отфильтрованные данные.
    """
    return df[df[column] >= min_phred]


def filter_by_read_length(df, column="read_length", min_length=50, max_length=None):
    """
    Удаляем риды, чья длина выходит за пределы [min_length, max_length].
    Если max_length=None, то фильтруем только по нижней границе.
    """
    res = df[df[column] >= min_length]
    if max_length is not None:
        res = res[res[column] <= max_length]
    return res


def filter_by_adapter(df, seq_column="sequence", adapters=None):
    """
    Удаляем риды, в которых обнаружены остатки адаптеров (последовательностей).
    Логика: если рид содержит любую подстроку из 'adapters' → рид выбрасываем.
    
    Параметры:
      seq_column (str): имя столбца, где хранится строка нуклеотидов.
      adapters (list): список строк, которые считаются адаптерами.
      
    Возвращает:
      DataFrame: риды, не содержащие ни одного адаптера.
    """
    if not adapters:
        # если адаптеры не заданы, просто возвращаем всё как есть
        return df

    mask = []
    for seq in df[seq_column]:
        # проверяем, содержит ли seq хотя бы один адаптер
        has_adapter = any(adpt in seq for adpt in adapters)
        mask.append(not has_adapter)  
    mask = pd.Series(mask, index=df.index)
    return df[mask]


def filter_by_n_count(df, seq_column="sequence", max_n=0):
    """
    Удаляем риды, где кол-во 'N' в последовательности > max_n.
    По умолчанию max_n=0 → удаляем все, у кого есть хотя бы один 'N'.
    """
    def has_too_many_n(seq):
        return seq.count("N") > max_n

    mask = df[seq_column].apply(lambda s: not has_too_many_n(s))
    return df[mask]


def filter_by_coverage(df, column="coverage", min_cov=10):
    """
    Удаляем строки (например, регионы, риды) с coverage < min_cov.
    """
    return df[df[column] >= min_cov]


def filter_by_mapping(df, mapq_col="mapq", mapq_threshold=20,
                      multi_mapping_col="num_mappings", remove_multimapped=True):
    """
    Удаляем риды, где MAPQ < mapq_threshold, 
    а при remove_multimapped=True – и те, где num_mappings > 1 (считаем многократно замапленные некачественными).
    """
    df_filtered = df[df[mapq_col] >= mapq_threshold]

    if remove_multimapped and multi_mapping_col in df_filtered.columns:
        df_filtered = df_filtered[df_filtered[multi_mapping_col] <= 1]
    return df_filtered


def remove_pcr_duplicates(df, chr_col="chr", start_col="start", end_col="end", strand_col="strand"):
    """
    Удаляем PCR-дубликаты по наивному критерию:
      (chr, start, end, strand) – уникальный "ключ". Все риды с одинаковым ключом, кроме первого, удаляем.
    """
    # строим tuple-ключ
    keys = df[[chr_col, start_col, end_col, strand_col]].astype(str).apply(tuple, axis=1)
    seen = set()
    to_keep = []

    for idx, key in zip(df.index, keys):
        if key not in seen:
            seen.add(key)
            to_keep.append(idx)

    return df.loc[to_keep]


def filter_by_gc_content(df, seq_column="sequence", min_gc=0.0, max_gc=1.0):
    """
    Удаляем риды, у которых GC-доля за пределами [min_gc, max_gc].
    Если min_gc, max_gc ≤ 1.0 => это доля (0..1). Иначе => проценты (0..100).
    """
    is_fraction = (max_gc <= 1.0 and min_gc >= 0.0)

    def gc_fraction(seq):
        gc_count = sum(base in "GCgc" for base in seq)
        return gc_count / len(seq) if len(seq) else 0

    def check_gc(seq):
        frac = gc_fraction(seq)
        if is_fraction:
            return (frac >= min_gc) and (frac <= max_gc)
        else:
            pct = frac * 100
            return (pct >= min_gc) and (pct <= max_gc)

    mask = df[seq_column].apply(check_gc)
    return df[mask]
