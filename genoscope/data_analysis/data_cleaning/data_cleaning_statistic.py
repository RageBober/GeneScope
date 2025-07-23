"""data_analysis.data_cleaning.data_cleaning_statistic
====================================================
Основные статистические приёмы очистки табличных данных:
• удаление дубликатов;
• заполнение NaN (ffill / bfill / mean / median / mode / interpolate / ml*);
• детекция выбросов (IQR, Z‑score, Mahalanobis).

*ML‑импутация делегируется модулю :pymod:`data_analysis.data_cleaning_ml` —
статистический код остаётся «лёгким» и не тянет `sklearn` при импорте.
"""

from __future__ import annotations

# ── STD / THIRD‑PARTY ───────────────────────────────────────────────
import logging
from typing import Literal, Sequence

import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.covariance import EmpiricalCovariance

# ── LOCAL ────────────────────────────────────────────────────────────
from genoscope.utils.columns import numeric  # helper → возвращает только числовые колонки

logger = logging.getLogger(__name__)
DF = pd.DataFrame

# ════════════════════════════════════════════════════════════════════
# 1. Дубликаты
# ════════════════════════════════════════════════════════════════════

def remove_duplicates(df: DF) -> DF:
    """Удаляет *полностью* одинаковые строки и логирует число удалённых."""
    if df.empty:
        logger.warning("[remove_duplicates] Пустой DataFrame — возврат без изменений")
        return df
    out = df.drop_duplicates()
    logger.info("[remove_duplicates] удалено %d дубликатов", len(df) - len(out))
    return out

# ════════════════════════════════════════════════════════════════════
# 2. Заполнение пропусков (NaN)
# ════════════════════════════════════════════════════════════════════

def handle_missing_values(
    df: DF,
    method: str = "ffill",
    columns: Sequence[str] | None = None,
    *,
    allow_mixed_types: bool = False,
) -> DF:
    if df.empty:
        logger.warning("[handle_missing_values] Пустой DataFrame — возврат без изменений")
        return df

    if method == "dialog":
        method = ask_user_method()
        logger.info("[handle_missing_values] выбран через dialog → %s", method)

    if columns is None:
        columns = list(df.columns)

    # игнорируем несуществующие колонки
    columns = [c for c in columns if c in df.columns]
    if not columns:
        return df

    dfc = df.copy()

    try:
        # ―― Простые методы pandas ――――――――――――――――――――――
        if method == "ffill":
            dfc[columns] = dfc[columns].ffill()
        elif method == "bfill":
            dfc[columns] = dfc[columns].bfill()

        # ―― Mean / Median ―――――――――――――――――――――――――――――
        elif method == "mean":
            num = numeric(dfc, columns, allow_mixed_types, ctx="mean")
            means = dfc[num].mean(numeric_only=True)
            dfc[num] = dfc[num].fillna(means)
        elif method == "median":
            
            num = numeric(dfc, columns, allow_mixed_types, ctx="median")
            med = dfc[num].median(numeric_only=True)
            dfc[num] = dfc[num].fillna(med)

        # ―― Mode ――――――――――――――――――――――――――――――――――――――
        elif method == "mode":
            for col in columns:
                m = dfc[col].mode(dropna=True)
                if not m.empty:
                    dfc[col] = dfc[col].fillna(m.iloc[0])

        # ―― Interpolate ――――――――――――――――――――――――――――――
        elif method == "interpolate":
            num = numeric(dfc, columns, allow_mixed_types, ctx="interpolate")
            dfc[num] = dfc[num].interpolate(method="linear")

        # ―― ML‑импутация ―――――――――――――――――――――――――――――――
        elif method == "ml":
            from data_cleaning.data_cleaning_ml import fill_missing_with_ml
            dfc = fill_missing_with_ml(dfc, columns)

        else:
            raise ValueError(method)

    except Exception as exc:  # noqa: BLE001
        logger.error("[handle_missing_values] %s", exc, exc_info=True)
        return df

    logger.info(
        "[handle_missing_values] '%s' завершён; осталось NaN: %d",
        method,
        int(dfc[columns].isna().sum().sum()),
    )
    return dfc

# ────────────────────────────────────────────────────────────────────

def ask_user_method() -> str:  # pragma: no cover – интерактив
    choices = ["mean", "median", "mode", "interpolate", "ffill", "bfill", "ml"]
    for i, m in enumerate(choices, 1):
        print(f"{i}. {m}")
    sel = input("№ метода (Enter = 1): ") or "1"
    try:
        idx = int(sel)
        return choices[idx - 1] if 1 <= idx <= len(choices) else "mean"
    except ValueError:
        return "mean"

# ════════════════════════════════════════════════════════════════════
# 3. Выбросы (IQR / Z‑score / Mahalanobis)
# ════════════════════════════════════════════════════════════════════

def detect_outliers(
    df: DF,
    column: str | None = None,
    *,
    method: Literal["iqr", "z-score", "mahalanobis"] = "iqr",
    threshold: float = 1.5,
) -> DF | pd.Series[bool]:
    if df.empty:
        empty = pd.DataFrame(False, index=df.index, columns=df.columns)
        return empty[column] if column else empty

    num = df.select_dtypes(include=[np.number]).columns
    mask_df = pd.DataFrame(False, index=df.index, columns=num)

    try:
        if method == "iqr":
            for col in ([column] if column else num):
                q1, q3 = df[col].quantile([0.25, 0.75])
                iqr = q3 - q1
                lb, ub = q1 - threshold * iqr, q3 + threshold * iqr
                mask_df[col] = (df[col] < lb) | (df[col] > ub)

        elif method == "z-score":
            for col in ([column] if column else num):
                z = zscore(df[col].astype(float), nan_policy="omit")
                mask_df[col] = np.abs(z) > threshold

        elif method == "mahalanobis":
            sub = df[num].dropna()
            if len(sub) >= 2:
                cov = EmpiricalCovariance().fit(sub)
                dist = cov.mahalanobis(sub)
                mask_df.loc[sub.index[dist > threshold], :] = True
        else:
            raise ValueError(method)

    except Exception as exc:  # noqa: BLE001
        logger.exception("[detect_outliers] %s", exc)

    if column:
        return mask_df[column]
    return mask_df if method in {"iqr", "z-score"} else mask_df.any(axis=1)

# ════════════════════════════════════════════════════════════════════
__all__: list[str] = [
    "remove_duplicates",
    "handle_missing_values",
    "detect_outliers",
]
