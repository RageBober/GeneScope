"""data_analysis.data_cleaning.data_cleaning_ml
--------------------------------------------
Функции для ML-обработки табличных данных: заполнение пропусков
(RandomForest, LogisticRegression), детекция выбросов (IsolationForest).
См. также: data_cleaning_statistic для статистических методов.
"""

from __future__ import annotations

import logging
from typing import Literal


from sklearn.ensemble import IsolationForest, RandomForestRegressor
from sklearn.linear_model import LogisticRegression

import pandas as pd
import numpy as np


logger = logging.getLogger(__name__)
# ──────────────────────────────────────────────────────────────
# Тип-алиасы
# ──────────────────────────────────────────────────────────────


def fill_missing_with_ml(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        nan_mask = out[col].isna()
        if not nan_mask.any():
            continue
        # mypy-friendly: второй аргумент `.loc` → list[str]
        other_cols: list[str] = list(out.columns.difference([col]))
        X_other = out.loc[~nan_mask, other_cols]
        y_known = out.loc[~nan_mask, col]
        if X_other.empty:
            continue
        is_num = pd.api.types.is_numeric_dtype(out[col])
        model = (
            RandomForestRegressor(n_estimators=20, random_state=42)
            if is_num
            else LogisticRegression(max_iter=1000)
        )
        try:
            model.fit(X_other, y_known)
            preds = model.predict(out.loc[nan_mask, other_cols]) 
            out.loc[nan_mask, col] = preds
        except Exception as exc:
            logger.warning("[fill_missing_with_ml] %s – оставил NaN", exc)
    return out

def detect_outliers(
    df: pd.DataFrame,
    column: str | None = None,
    *,
    method: Literal["isolation_forest"] = "isolation_forest",
    threshold: float = 0.05,
) -> pd.DataFrame | pd.Series[bool]:
    if df.empty:
        empty = pd.DataFrame(False, index=df.index, columns=df.columns)
        return empty[column] if column else empty
    num_cols = df.select_dtypes(include=[np.number]).columns
    mask_df = pd.DataFrame(False, index=df.index, columns=num_cols)  
    try:
        if method == "isolation_forest":
            num = df[num_cols].dropna()
            if not num.empty:
                cont = threshold if 0 < threshold < 1 else 0.05
                iso = IsolationForest(contamination=cont, random_state=42)
                preds = iso.fit_predict(num)
                bad_idx = num.index[preds == -1]
                mask_df.loc[bad_idx, :] = True 
        else:
            raise ValueError(f"Unsupported method '{method}'")
    except Exception as exc:
        logger.exception("[detect_outliers] непредвидённая ошибка: %s", exc)
    if column:
        return mask_df[column]
    return mask_df.any(axis=1)