from __future__ import annotations

import traceback
from typing import Any, Dict, Literal, Sequence

import pandas as pd
from genoscope.data_analysis._typing import BoolS

# ---------------------------------------------------------------------------#
#  Тип-алиасы (можно вынести в отдельный types.py, если понадобится)
# ---------------------------------------------------------------------------#
DF = pd.DataFrame

PipelineStep = Dict[str, Any]
Pipeline = Sequence[PipelineStep]

##############################################################################
#  1.  Базовые «атомарные» фильтры
##############################################################################


def _col_missing(df: DF, col: str, ctx: str) -> bool:
    """Возвращает True, если столбца нет (и печатает предупреждение)."""
    if col not in df.columns:
        print(f"[{ctx}] Колонка «{col}» не найдена, пропускаем фильтр.")
        return True
    return False


def filter_by_qual(df: DF, *, qual_col: str = "QUAL", min_qual: float = 30.0) -> DF:
    """Оставляет варианты с QUAL ≥ min_qual."""
    if _col_missing(df, qual_col, "filter_by_qual"):
        return df
    mask: BoolS = df[qual_col].astype(float) >= min_qual
    return df.loc[mask]


def filter_by_depth(df: DF, *, depth_col: str = "DP", min_depth: int = 10) -> DF:
    """Оставляет варианты с глубиной покрытия DP ≥ min_depth."""
    if _col_missing(df, depth_col, "filter_by_depth"):
        return df
    mask: BoolS = df[depth_col].astype(float) >= min_depth
    return df.loc[mask]


def filter_by_strand_bias(
    df: DF,
    *,
    sb_col: str = "SB",
    max_bias: float = 0.01,
) -> DF:
    """
    Пример фильтра по strand-bias: оставляем строки, где SB > max_bias
    (подразумеваем: «чем меньше число, тем сильнее артефакт»).
    """
    if _col_missing(df, sb_col, "filter_by_strand_bias"):
        return df
    mask: BoolS = df[sb_col].astype(float) > max_bias
    return df.loc[mask]


def filter_by_allele_frequency(
    df: DF,
    *,
    af_col: str = "AF",
    min_af: float = 0.01,
    max_af: float = 0.99,
) -> DF:
    """Оставляет варианты c min_af ≤ AF ≤ max_af."""
    if _col_missing(df, af_col, "filter_by_allele_frequency"):
        return df
    af = df[af_col].astype(float)
    mask: BoolS = (af >= min_af) & (af <= max_af)
    return df.loc[mask]


def filter_by_annotation(
    df: DF,
    *,
    ann_col: str = "ANN",
    allowed_types: Sequence[str] | None = None,
) -> DF:
    """
    Удаляет варианты, НЕ содержащие ни одну из строк из `allowed_types`
    в поле `ANN`.
    """
    if not allowed_types:
        return df
    if _col_missing(df, ann_col, "filter_by_annotation"):
        return df

    def has_allowed(row_val: str) -> bool:
        rv = str(row_val)
        return any(a_type in rv for a_type in allowed_types)

    mask: BoolS = df[ann_col].apply(has_allowed)
    return df.loc[mask]


def filter_by_db_snp_clinvar(
    df: DF,
    *,
    db_col: str = "dbSNP_ID",
    clinvar_col: str = "ClinVar_ID",
    require_dbsnp: bool = True,
    require_clinvar: bool = False,
) -> DF:
    """Фильтрация по наличию ID в dbSNP / ClinVar."""
    res = df
    if require_dbsnp:
        if _col_missing(res, db_col, "filter_by_db_snp_clinvar-dbSNP"):
            pass
        else:
            res = res[res[db_col].notna()]

    if require_clinvar:
        if _col_missing(res, clinvar_col, "filter_by_db_snp_clinvar-ClinVar"):
            pass
        else:
            res = res[res[clinvar_col].notna()]

    return res


##############################################################################
#  2.  Пайплайн «variant-calling» фильтров
##############################################################################

StepType = Literal[
    "qual",
    "depth",
    "strand_bias",
    "allele_freq",
    "annotation",
    "db_snp_clinvar",
]


def apply_variant_filter_pipeline(df: DF, pipeline: Pipeline) -> DF:
    """
    Применяет список шагов вида::

        {"type": "qual", "qual_col": "QUAL", "min_qual": 30.0}

    Неизвестный `type` → печать предупреждения и пропуск шага.
    Ошибка внутри шага → лог + продолжаем пайплайн (не прерываем).
    """
    for step in pipeline:
        step_type: StepType | str | None = step.get("type")
        try:
            if step_type == "qual":
                df = filter_by_qual(
                    df,
                    qual_col=step.get("qual_col", "QUAL"),
                    min_qual=step.get("min_qual", 30.0),
                )

            elif step_type == "depth":
                df = filter_by_depth(
                    df,
                    depth_col=step.get("depth_col", "DP"),
                    min_depth=step.get("min_depth", 10),
                )

            elif step_type == "strand_bias":
                df = filter_by_strand_bias(
                    df,
                    sb_col=step.get("sb_col", "SB"),
                    max_bias=step.get("max_bias", 0.01),
                )

            elif step_type == "allele_freq":
                df = filter_by_allele_frequency(
                    df,
                    af_col=step.get("af_col", "AF"),
                    min_af=step.get("min_af", 0.01),
                    max_af=step.get("max_af", 0.99),
                )

            elif step_type == "annotation":
                df = filter_by_annotation(
                    df,
                    ann_col=step.get("ann_col", "ANN"),
                    allowed_types=step.get("allowed_types"),
                )

            elif step_type == "db_snp_clinvar":
                df = filter_by_db_snp_clinvar(
                    df,
                    db_col=step.get("db_col", "dbSNP_ID"),
                    clinvar_col=step.get("clinvar_col", "ClinVar_ID"),
                    require_dbsnp=step.get("require_dbsnp", True),
                    require_clinvar=step.get("require_clinvar", False),
                )

            else:
                print(f"[apply_variant_filter_pipeline] Неизвестный type: {step_type}")

        except Exception as exc:  # noqa: BLE001
            print(f"[apply_variant_filter_pipeline] Ошибка на шаге «{step_type}»: {exc}")
            traceback.print_exc()

    return df


##############################################################################
#  3.  Черновой «однокнопочный» фильтр
##############################################################################


def filter_variants(
    df: DF,
    *,
    qual_col: str = "QUAL",
    min_qual: float = 30.0,
) -> DF:
    """Простейший фильтр: QUAL ≥ min_qual."""
    if _col_missing(df, qual_col, "filter_variants"):
        return df
    return df.loc[df[qual_col].astype(float) >= min_qual]


__all__: Sequence[str] = [
    "filter_by_qual",
    "filter_by_depth",
    "filter_by_strand_bias",
    "filter_by_allele_frequency",
    "filter_by_annotation",
    "filter_by_db_snp_clinvar",
    "apply_variant_filter_pipeline",
    "filter_variants",
]
