import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.covariance import EmpiricalCovariance
from sklearn.ensemble import IsolationForest


def get_outlier_mask(
    df: pd.DataFrame,
    column: str | None = None,
    method: str = "iqr",
    threshold: float = 1.5,
) -> pd.Series | pd.DataFrame:
    """
    Возвращает булеву маску выбросов в DataFrame по выбранному методу.

    Параметры:
        df (pd.DataFrame): Исходные данные.
        column (str или None): Имя столбца для анализа (или None для всех числовых столбцов).
        method (str): Метод обнаружения выбросов ('iqr', 'z-score', 'mahalanobis', 'isolation_forest').
        threshold (float): Пороговый параметр для выбранного метода.

    Возвращает:
        pd.Series или pd.DataFrame: Булева маска (True — выброс, False — нормальное значение).
    """
    if df.empty:
        return pd.DataFrame(False, index=df.index, columns=[column] if column else df.columns)

    if method in ["iqr", "z-score"]:
        cols = [column] if column else df.select_dtypes(include=[np.number]).columns
        mask = pd.DataFrame(False, index=df.index, columns=cols)
        for col in cols:
            if method == "iqr":
                q1 = df[col].quantile(0.25)
                q3 = df[col].quantile(0.75)
                iqr = q3 - q1
                lb = q1 - threshold * iqr
                ub = q3 + threshold * iqr
                mask[col] = (df[col] < lb) | (df[col] > ub)
            elif method == "z-score":
                z_scores = zscore(df[col].astype(float), nan_policy="omit")
                mask[col] = np.abs(z_scores) > threshold
        return mask[column] if column else mask

    elif method == "mahalanobis":
        # Для mahalanobis только по всем числовым столбцам!
        num_df = df.select_dtypes(include=[np.number]).dropna(axis=0)
        if num_df.shape[1] < 2:
            # Махаланобис не применим к одному столбцу
            return pd.Series(False, index=df.index)
        cov = EmpiricalCovariance().fit(num_df)
        distances = cov.mahalanobis(num_df)
        msk = distances > threshold
        mask = pd.Series(False, index=df.index)
        mask.loc[num_df.index] = msk
        return mask

    elif method == "isolation_forest":
        # Скоро перенесём в ML, пока оставляем здесь для совместимости
        cont = float(threshold) if 0 < threshold < 1 else 0.05
        num_df = df.select_dtypes(include=[np.number]).dropna(axis=0)
        if num_df.shape[0] < 10:
            # Нет смысла применять лес на слишком малом количестве данных
            return pd.Series(False, index=df.index)
        iso = IsolationForest(contamination=cont, random_state=0)
        preds = iso.fit_predict(num_df)  # +1 норм, -1 выброс
        msk = preds == -1
        mask = pd.Series(False, index=df.index)
        mask.loc[num_df.index] = msk
        return mask

    else:
        raise ValueError(f"Unknown method: {method}")
