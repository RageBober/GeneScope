from __future__ import annotations

import traceback
from collections.abc import Sequence
from typing import Any
from typing import Literal


from genoscope.data_analysis._typing import BoolS, DF

# ──────────────────────────────────────────────────────────────────────────────
# Type aliases
# ──────────────────────────────────────────────────────────────────────────────
PipelineStep = dict[str, Any]
Pipeline = Sequence[PipelineStep]


# ──────────────────────────────────────────────────────────────────────────────
# 1. Атомарные фильтры
# ──────────────────────────────────────────────────────────────────────────────
def _col_missing(df: DF, col: str, ctx: str) -> bool:
    if col not in df.columns:
        print(f"[{ctx}] Колонка «{col}» не найдена, пропускаем.")
        return True
    return False


def filter_by_otu_coverage(
    df: DF,
    *,
    coverage_col: str = "OTU_COUNT",
    min_coverage: int = 100,
) -> DF:
    if _col_missing(df, coverage_col, "filter_by_otu_coverage"):
        return df
    mask: BoolS = df[coverage_col].astype(float) >= min_coverage
    return df.loc[mask]


def filter_rare_taxa(
    df: DF,
    *,
    taxon_col: str = "TAXON",
    sample_col: str = "SAMPLE_COUNT",
    min_samples: int = 2,
) -> DF:
    if _col_missing(df, sample_col, "filter_rare_taxa"):
        return df
    mask: BoolS = df[sample_col].astype(int) >= min_samples
    return df.loc[mask]


def filter_by_maf_epigenetics(
    df: DF,
    *,
    maf_col: str = "MAF",
    min_maf: float = 0.05,
) -> DF:
    if _col_missing(df, maf_col, "filter_by_maf_epigenetics"):
        return df
    mask: BoolS = df[maf_col].astype(float) >= min_maf
    return df.loc[mask]


def filter_by_cpg_confidence(
    df: DF,
    *,
    coverage_col: str = "CpG_COV",
    pvalue_col: str = "CpG_PVALUE",
    min_coverage: int = 10,
    max_pvalue: float = 0.05,
) -> DF:
    res = df
    if not _col_missing(res, coverage_col, "filter_by_cpg_confidence-coverage"):
        res = res.loc[res[coverage_col].astype(float) >= min_coverage]

    if not _col_missing(res, pvalue_col, "filter_by_cpg_confidence-pvalue"):
        res = res.loc[res[pvalue_col].astype(float) <= max_pvalue]

    return res


# ──────────────────────────────────────────────────────────────────────────────
# 2. Пайплайн
# ──────────────────────────────────────────────────────────────────────────────
StepType = Literal[
    "otu_coverage",
    "rare_taxa",
    "maf_epigenetics",
    "cpg_confidence",
]


def apply_meta_epi_filter_pipeline(
    df: DF,
    pipeline: Pipeline,
) -> DF:
    """
    Применяет последовательность шагов к таблице мета/эпигенетики.
    Неизвестный type → пропуск шага. Ошибка в шаге → сообщение и продолжение.
    """
    for step in pipeline:
        step_type: StepType | str | None = step.get("type")
        try:
            if step_type == "otu_coverage":
                df = filter_by_otu_coverage(
                    df,
                    coverage_col=step.get("coverage_col", "OTU_COUNT"),
                    min_coverage=step.get("min_coverage", 100),
                )

            elif step_type == "rare_taxa":
                df = filter_rare_taxa(
                    df,
                    taxon_col=step.get("taxon_col", "TAXON"),
                    sample_col=step.get("sample_col", "SAMPLE_COUNT"),
                    min_samples=step.get("min_samples", 2),
                )

            elif step_type == "maf_epigenetics":
                df = filter_by_maf_epigenetics(
                    df,
                    maf_col=step.get("maf_col", "MAF"),
                    min_maf=step.get("min_maf", 0.05),
                )

            elif step_type == "cpg_confidence":
                df = filter_by_cpg_confidence(
                    df,
                    coverage_col=step.get("coverage_col", "CpG_COV"),
                    pvalue_col=step.get("pvalue_col", "CpG_PVALUE"),
                    min_coverage=step.get("min_coverage", 10),
                    max_pvalue=step.get("max_pvalue", 0.05),
                )

            else:
                print(f"[apply_meta_epi_filter_pipeline] Неизвестный type: {step_type}")

        except Exception as exc:
            print(
                f"[apply_meta_epi_filter_pipeline] Ошибка на шаге «{step_type}»: {exc}"
            )
            traceback.print_exc()

    return df


# ──────────────────────────────────────────────────────────────────────────────
__all__: Sequence[str] = [
    "apply_meta_epi_filter_pipeline",
    "filter_by_cpg_confidence",
    "filter_by_maf_epigenetics",
    "filter_by_otu_coverage",
    "filter_rare_taxa",
]
