# utils/outlier_utils.py
from __future__ import annotations

from typing import Literal
from typing import cast
from typing import overload

import numpy as np
import pandas as pd
from numpy.typing import NDArray
from scipy.stats import zscore
from sklearn.covariance import EmpiricalCovariance
from sklearn.ensemble import IsolationForest

BoolSeries = pd.Series[bool]  # алиас удобнее читать
BoolFrame = pd.DataFrame  # pandas-stubs ещё не умеет DataFrame[bool]

MethodT = Literal["iqr", "z-score", "mahalanobis", "isolation_forest"]


@overload
def get_outlier_mask(
    df: pd.DataFrame,
    *,
    column: str,
    method: MethodT = "iqr",
    threshold: float = 1.5,
) -> BoolSeries: ...  # одиночный столбец → Series[bool]


@overload
def get_outlier_mask(
    df: pd.DataFrame,
    *,
    column: None = None,
    method: MethodT = "iqr",
    threshold: float = 1.5,
) -> BoolFrame: ...  # все числовые столбцы → DataFrame


def get_outlier_mask(
    df: pd.DataFrame,
    *,
    column: str | None = None,
    method: MethodT = "iqr",
    threshold: float = 1.5,
) -> BoolSeries | BoolFrame:
    """
    Возвращает булеву маску выбросов.

    Args:
        df: исходный DataFrame.
        column: имя столбца или None (по всем числовым).
        method: 'iqr' | 'z-score' | 'mahalanobis' | 'isolation_forest'.
        threshold: параметр порога для выбранного метода.

    Returns:
        Series[bool] — если указан column,
        DataFrame     — если column is None.
    """
    # --- защита от пустого входа ------------------------------------------
    if df.empty:
        empty_mask: BoolFrame = pd.DataFrame(
            False,
            index=df.index,
            columns=[column] if column else df.columns,
        )
        return cast(
            "BoolSeries | BoolFrame",  # mypy OK
            empty_mask[column] if column else empty_mask,
        )

    # --- IQR / Z-score -----------------------------------------------------
    if method in ("iqr", "z-score"):
        cols = [column] if column else list(df.select_dtypes(np.number).columns)
        mask_df: BoolFrame = pd.DataFrame(False, index=df.index, columns=cols)

        for col in cols:
            if method == "iqr":
                q1, q3 = df[col].quantile([0.25, 0.75])
                iqr = q3 - q1
                lb, ub = q1 - threshold * iqr, q3 + threshold * iqr
                mask_df[col] = (df[col] < lb) | (df[col] > ub)

            else:  # 'z-score'
                z: NDArray[np.float_] = zscore(df[col].astype(float), nan_policy="omit")
                mask_df[col] = np.abs(z) > threshold

        return cast("BoolSeries | BoolFrame", mask_df[cols[0]] if column else mask_df)

    # --- Mahalanobis -------------------------------------------------------
    if method == "mahalanobis":
        num_df = df.select_dtypes(np.number).dropna()

        # метод не применим к 0/1 столбцу
        if num_df.shape[1] < 2:
            return pd.Series(False, index=df.index)

        cov = EmpiricalCovariance().fit(num_df)
        distances = cov.mahalanobis(num_df)
        flag: NDArray[np.bool_] = distances > threshold

        mask = pd.Series(False, index=df.index)
        mask.loc[num_df.index] = flag
        return mask

    # --- Isolation Forest --------------------------------------------------
    if method == "isolation_forest":
        cont = threshold if 0 < threshold < 1 else 0.05
        num_df = df.select_dtypes(np.number).dropna()

        if num_df.shape[0] < 10:
            return pd.Series(False, index=df.index)

        iso = IsolationForest(contamination=cont, random_state=0)
        preds: NDArray[np.int_] = iso.fit_predict(num_df)  # + 1 норм, -1 выброс
        flag = preds == -1

        mask = pd.Series(False, index=df.index)
        mask.loc[num_df.index] = flag
        return mask

    # --- Неподдерживаемый метод -------------------------------------------
    raise ValueError(f"Unknown method: {method!r}")
