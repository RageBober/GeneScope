from __future__ import annotations

import traceback
from collections.abc import Sequence
from typing import Any
from typing import Literal
from typing import cast

import numpy as np
import pandas as pd

from genoscope.data_analysis._typing import BoolS, DF

# ──────────────────────────────────────────────────────────────────────────────
# type aliases
# ──────────────────────────────────────────────────────────────────────────────


PipeStep = dict[str, Any]
Pipeline = Sequence[PipeStep]


# ──────────────────────────────────────────────────────────────────────────────
# helpers
# ──────────────────────────────────────────────────────────────────────────────
def _col_missing(df: DF, col: str, ctx: str) -> bool:
    if col not in df.columns:
        print(f"[{ctx}] Колонка «{col}» не найдена, пропускаем.")
        return True
    return False


# ──────────────────────────────────────────────────────────────────────────────
# 1. атомарные фильтры
# ──────────────────────────────────────────────────────────────────────────────
def filter_by_expression_level(
    df: DF,
    *,
    column: str = "TPM",
    threshold: float = 1.0,
    how: Literal["above", "below"] = "above",
) -> DF:
    if _col_missing(df, column, "filter_by_expression_level"):
        return df

    series = df[column].astype(float)
    mask: BoolS = series >= threshold if how == "above" else series <= threshold
    return df.loc[mask]


def filter_by_zero_count(
    df: DF,
    *,
    columns: Sequence[str] | None = None,
    max_zero_fraction: float = 0.5,
) -> DF:
    if columns is None:
        cols: Sequence[str] = list(df.select_dtypes(include=[np.number]).columns)
    else:
        cols = tuple(columns)  # превращаем в immutable-Sequence для mypy

    if not cols:
        return df

    zero_frac: pd.Series[float] = (df[cols] == 0).sum(axis=1) / float(len(cols))
    mask: BoolS = zero_frac <= max_zero_fraction
    return df.loc[mask]


def filter_by_fold_change(
    df: DF,
    *,
    column_fc: str = "log2FC",
    min_fc: float = 1.0,
) -> DF:
    if _col_missing(df, column_fc, "filter_by_fold_change"):
        return df

    # было: np.abs(df[column_fc].astype(float))  -> ndarray
    series_fc = df[column_fc].astype(float).abs()  # ← Series
    mask: BoolS = series_fc >= min_fc
    return df.loc[mask]


def filter_by_significance(
    df: DF,
    *,
    column_p: str = "pvalue",
    alpha: float = 0.05,
    method: Literal["pvalue", "fdr"] = "pvalue",
) -> DF:
    if _col_missing(df, column_p, "filter_by_significance"):
        return df
    mask: BoolS = df[column_p].astype(float) <= alpha
    return df.loc[mask]


def filter_by_biomarker_bases(
    df: DF,
    *,
    gene_column: str = "gene_id",
    biomarker_list: Sequence[str] | None = None,
) -> DF:
    if biomarker_list is None:
        return df
    if _col_missing(df, gene_column, "filter_by_biomarker_bases"):
        return df
    mask: BoolS = df[gene_column].isin(biomarker_list)
    return df.loc[mask]


# ──────────────────────────────────────────────────────────────────────────────
# 2. apply-pipeline
# ──────────────────────────────────────────────────────────────────────────────
StepType = Literal[
    "expr_level",
    "zero_count",
    "fold_change",
    "significance",
    "biomarker",
]


def apply_expression_filters(df: DF, pipeline: Pipeline | None = None) -> DF:
    if pipeline is None:
        return df

    for step in pipeline:
        step_type: StepType | str | None = step.get("type")
        try:
            if step_type == "expr_level":
                df = filter_by_expression_level(
                    df,
                    column=step.get("column", "TPM"),
                    threshold=cast("float", step.get("threshold", 1.0)),
                    how=cast("Literal['above', 'below']", step.get("how", "above")),
                )

            elif step_type == "zero_count":
                df = filter_by_zero_count(
                    df,
                    columns=step.get("columns"),
                    max_zero_fraction=cast("float", step.get("max_zero_fraction", 0.5)),
                )

            elif step_type == "fold_change":
                df = filter_by_fold_change(
                    df,
                    column_fc=step.get("column_fc", "log2FC"),
                    min_fc=cast("float", step.get("min_fc", 1.0)),
                )

            elif step_type == "significance":
                df = filter_by_significance(
                    df,
                    column_p=step.get("column_p", "pvalue"),
                    alpha=cast("float", step.get("alpha", 0.05)),
                    method=cast(
                        "Literal['pvalue', 'fdr']", step.get("method", "pvalue")
                    ),
                )

            elif step_type == "biomarker":
                df = filter_by_biomarker_bases(
                    df,
                    gene_column=step.get("gene_column", "gene_id"),
                    biomarker_list=step.get("biomarker_list"),
                )

            else:
                print(f"[apply_expression_filters] Неизвестный type: {step_type}")

        except Exception as exc:
            print(f"[apply_expression_filters] Ошибка на шаге «{step_type}»: {exc}")
            traceback.print_exc()

    return df


# ──────────────────────────────────────────────────────────────────────────────
__all__: Sequence[str] = [
    "apply_expression_filters",
    "filter_by_biomarker_bases",
    "filter_by_expression_level",
    "filter_by_fold_change",
    "filter_by_significance",
    "filter_by_zero_count",
]
