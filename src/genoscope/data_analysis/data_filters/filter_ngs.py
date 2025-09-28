from __future__ import annotations

from collections.abc import Sequence
from typing import Any
from typing import cast

import pandas as pd

from genoscope.data_analysis._typing import BoolS, DF

# ──────────────────────────────────────────────────────────────────────────────
#  Тип-алиасы (опционально можно вынести в общий types.py)
# ──────────────────────────────────────────────────────────────────────────────


# ──────────────────────────────────────────────────────────────────────────────
# 1. Качество и длина ридов
# ──────────────────────────────────────────────────────────────────────────────
def filter_by_phred_quality(
    df: DF,
    *,
    column: str = "mean_phred",
    min_phred: float = 20,
) -> DF:
    if column not in df.columns:
        return df
    mask: BoolS = df[column].astype(float) >= min_phred
    return df.loc[mask]


def filter_by_read_length(
    df: DF,
    *,
    column: str = "read_length",
    min_length: int = 50,
    max_length: int | None = None,
) -> DF:
    if column not in df.columns:
        return df
    res = df.loc[df[column].astype(int) >= min_length]
    if max_length is not None:
        res = res.loc[res[column].astype(int) <= max_length]
    return res


# ──────────────────────────────────────────────────────────────────────────────
# 2. Адаптеры, “N”-символы, coverage
# ──────────────────────────────────────────────────────────────────────────────
def filter_by_adapter(
    df: DF,
    *,
    seq_column: str = "sequence",
    adapters: Sequence[str] | None = None,
) -> DF:
    if not adapters or seq_column not in df.columns:
        return df
    mask: BoolS = (
        df[seq_column]
        .astype(str)
        .apply(lambda s: not any(adpt in s for adpt in adapters))
    )
    return df.loc[mask]


def filter_by_n_count(
    df: DF,
    *,
    seq_column: str = "sequence",
    max_n: int = 0,
) -> DF:
    if seq_column not in df.columns:
        return df

    def has_too_many_n(seq: str) -> bool:
        return seq.upper().count("N") > max_n

    mask: BoolS = df[seq_column].astype(str).apply(lambda s: not has_too_many_n(s))
    return df.loc[mask]


def filter_by_coverage(
    df: DF,
    *,
    column: str = "coverage",
    min_cov: float = 10,
) -> DF:
    if column not in df.columns:
        return df
    return df.loc[df[column].astype(float) >= min_cov]


# ──────────────────────────────────────────────────────────────────────────────
# 3. MAPQ, multi-mapping
# ──────────────────────────────────────────────────────────────────────────────
def filter_by_mapping(
    df: DF,
    *,
    mapq_col: str = "mapq",
    mapq_threshold: int = 20,
    multi_mapping_col: str = "num_mappings",
    remove_multimapped: bool = True,
) -> DF:
    if mapq_col not in df.columns:
        return df

    df_filtered = df.loc[df[mapq_col].astype(int) >= mapq_threshold]

    if remove_multimapped and multi_mapping_col in df_filtered.columns:
        df_filtered = df_filtered.loc[df_filtered[multi_mapping_col].astype(int) <= 1]

    return df_filtered


# ──────────────────────────────────────────────────────────────────────────────
# 4. PCR-дубликаты
# ──────────────────────────────────────────────────────────────────────────────
def remove_pcr_duplicates(
    df: DF,
    *,
    chr_col: str = "chr",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: str = "strand",
) -> DF:
    required = {chr_col, start_col, end_col, strand_col}
    if not required.issubset(df.columns):
        return df

    keys: pd.Series[Any] = (
        df[[chr_col, start_col, end_col, strand_col]].astype(str).apply(tuple, axis=1)
    )
    seen: set[tuple[str, str, str, str]] = set()
    keep_idx: list[int] = []

    for idx, key_obj in zip(df.index.to_list(), keys):

        key = cast("tuple[str, str, str, str]", key_obj)
        if key not in seen:
            seen.add(key)
            keep_idx.append(idx)

    return df.loc[keep_idx]


# ──────────────────────────────────────────────────────────────────────────────
# 5. GC-контент
# ──────────────────────────────────────────────────────────────────────────────
def _gc_fraction(seq: str) -> float:
    seq_up = seq.upper()
    gc = sum(base in "GC" for base in seq_up)
    return gc / len(seq_up) if seq_up else 0.0


def filter_by_gc_content(
    df: DF,
    *,
    seq_column: str = "sequence",
    min_gc: float = 0.0,
    max_gc: float = 1.0,
) -> DF:
    if seq_column not in df.columns:
        return df

    is_fraction = 0.0 <= min_gc <= 1.0 and 0.0 <= max_gc <= 1.0

    def inside_bounds(seq: str) -> bool:
        frac = _gc_fraction(seq)
        val = frac if is_fraction else frac * 100
        return min_gc <= val <= max_gc

    mask: BoolS = df[seq_column].astype(str).apply(inside_bounds)
    return df.loc[mask]


# ──────────────────────────────────────────────────────────────────────────────
# 6. Финальный «в лоб» фильтр числа ридов
# ──────────────────────────────────────────────────────────────────────────────
def filter_ngs(
    df: DF,
    *,
    min_reads: int = 10,
    column: str = "READS",
) -> DF:
    """Оставляет строки, где `column` ≥ min_reads."""
    if column not in df.columns:
        print(f"[filter_ngs] Колонка «{column}» не найдена — пропустил фильтр.")
        return df
    return df.loc[df[column].astype(int) >= min_reads]


# ──────────────────────────────────────────────────────────────────────────────
__all__: Sequence[str] = [
    "filter_by_adapter",
    "filter_by_coverage",
    "filter_by_gc_content",
    "filter_by_mapping",
    "filter_by_n_count",
    "filter_by_phred_quality",
    "filter_by_read_length",
    "filter_ngs",
    "remove_pcr_duplicates",
]
