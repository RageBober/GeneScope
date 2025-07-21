from __future__ import annotations

# ─── STDLIB ──────────────────────────────────────────────────────────
import logging
from typing import Literal

# ─── THIRD-PARTY ─────────────────────────────────────────────────────
import numpy as np
import pandas as pd
from scipy.stats import zscore
from sklearn.covariance import EmpiricalCovariance
from sklearn.ensemble import IsolationForest, RandomForestRegressor
from sklearn.linear_model import LogisticRegression

# ─── LOCAL / INTERNAL ───────────────────────────────────────────────
from genoscope.mlmodel.data_cleaning_AI.ml_imputation import fill_missing_with_ml

"""data_analysis.data_cleaning
---------------------------------
Функции для базовой пред‑ и пост‑обработки табличных (и частично —
последовательностных) данных: удаление дубликатов, заполнение NaN,
импутация ML, детекция выбросов. Файл содержит только «core‑MVP» —
никакой прод‑обвязки (Docker, Prometheus) и т. д.
"""


# ─── MODULE-LEVEL CONSTANTS & LOGGER ────────────────────────────────
logger = logging.getLogger(__name__)


def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """Удаляет дубликаты *полностью одинаковых* строк.

    Логирует количество удалённых строк, но всегда возвращает копию,
    чтобы не удивлять вызвавшего inplace‑модификациями.
    """
    if df.empty:
        logger.warning("[remove_duplicates] Пустой DataFrame — возврат без изменений")
        return df

    before = len(df)
    out = df.drop_duplicates()
    removed = before - len(out)
    if removed:
        logger.info("[remove_duplicates] удалено %d дубликатов", removed)
    return out


def handle_missing_values(
    df: pd.DataFrame,
    method: str = "ffill",
    columns: list = None,
    *,
    allow_mixed_types: bool = False,
) -> pd.DataFrame:
    """
    Обрабатывает пропущенные значения выбранными методами:
      'ffill', 'bfill', 'mean', 'median', 'mode', 'interpolate', 'ml',
      или 'dialog' — если хотим показать GUI для выбора метода.

    Параметры
    ---------
    df : pd.DataFrame
        Исходный DataFrame.
    method : str
        Метод, один из:
            'ffill', 'bfill', 'mean', 'median', 'mode',
            'interpolate', 'ml', 'dialog'.
    columns : list, optional
        Список столбцов, где заполнять пропуски; если None, все столбцы.
    allow_mixed_types : bool
        Если False, при 'mean'/'median'/'interpolate' нечисловые столбцы будут игнорироваться.

    Возвращает
    ----------
    pd.DataFrame
        Новый (скопированный) DataFrame с заполненными пропусками.
    """
    if df.empty:
        logger.warning("[handle_missing_values] Пустой DataFrame, возврат без изменений.")
        return df

    if method == "dialog":
        # Пример псевдо-GUI, где пользователь выбирает метод в консоли
        chosen = ask_user_method()  # Эта функция ниже — просто пример
        logger.info(f"[handle_missing_values] Пользователь выбрал метод через dialog: {chosen}")
        method = chosen

    if columns is None:
        columns = df.columns.tolist()

    # Проверим, что все указанные столбцы существуют
    missing_cols = [c for c in columns if c not in df.columns]
    if missing_cols:
        logger.warning(f"[handle_missing_values] Столбцы {missing_cols} не найдены. Пропускаем их.")
        columns = [c for c in columns if c in df.columns]
        if not columns:
            logger.warning(
                "[handle_missing_values] Ни одного подходящего столбца не осталось, возврат DF без изменений."
            )
            return df

    df_copy = df.copy()

    try:
        if method in ["ffill", "bfill"]:
            # Используем встроенные методы pandas fillna
            df_copy[columns] = df_copy[columns].fillna(method=method)

        elif method == "mean":
            # Для среднего берем только числовые столбцы
            num_cols = _check_numeric_columns(df_copy, columns, allow_mixed_types, context="mean")
            mean_vals = df_copy[num_cols].mean(numeric_only=True)
            df_copy[num_cols] = df_copy[num_cols].fillna(mean_vals)

        elif method == "median":
            # Медиана — только числовые столбцы
            num_cols = _check_numeric_columns(df_copy, columns, allow_mixed_types, context="median")
            median_vals = df_copy[num_cols].median(numeric_only=True)
            df_copy[num_cols] = df_copy[num_cols].fillna(median_vals)

        elif method == "mode":
            # Моду вычисляем для каждого столбца отдельно
            # (для числовых и категориальных без ограничений).
            for col in columns:
                series_mode = df_copy[col].mode(dropna=True)
                if not series_mode.empty:
                    # Если мод несколько, берем первый
                    fill_val = series_mode.iloc[0]
                    df_copy[col] = df_copy[col].fillna(fill_val)

        elif method == "interpolate":
            # Линейная интерполяция только для числовых
            num_cols = _check_numeric_columns(
                df_copy, columns, allow_mixed_types, context="interpolate"
            )
            df_copy[num_cols] = df_copy[num_cols].interpolate(method="linear")

        elif method == "ml":
            df_copy = fill_missing_with_ml(df_copy, columns)

        else:
            raise ValueError(f"[handle_missing_values] Unsupported method: '{method}'")

    except Exception as exc:
        logger.error(
            f"[handle_missing_values] Ошибка при заполнении методом '{method}': {exc}",
            exc_info=True,
        )
        return df  # Возвращаем исходный DataFrame при ошибке

    # Логируем, сколько пропусков осталось
    leftover_nans = df_copy[columns].isna().sum().sum()
    logger.info(
        f"[handle_missing_values] Метод '{method}' завершён. Осталось пропусков: {leftover_nans}"
    )
    return df_copy


# -------------------------------------------------------------------------
# Вспомогательная функция для проверки "числовых" столбцов
# -------------------------------------------------------------------------
def _check_numeric_columns(
    df: pd.DataFrame, columns: list, allow_mixed: bool, context: str
) -> list:
    numeric_cols = []
    for col in columns:
        if col not in df.columns:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            numeric_cols.append(col)
        else:
            if allow_mixed:
                logger.info(
                    f"[_check_numeric_columns] '{col}' нечисловой, но allow_mixed_types=True => пропускаем предупреждение."
                )
            else:
                logger.warning(
                    f"[_check_numeric_columns] '{col}' нечисловой, исключаем для метода '{context}'."
                )
    return numeric_cols


# -------------------------------------------------------------------------
# Пример ПСЕВДО-функции, как будто открываем GUI
# -------------------------------------------------------------------------
def ask_user_method() -> str:
    """
    Имитирует GUI-спрос у пользователя (в реальном коде это мог бы быть диалог Tkinter).
    Возвращает строку: 'mean' / 'median' / 'mode' / 'ml' / ...
    """

    possible_methods = ["mean", "median", "mode", "ml", "interpolate", "ffill", "bfill"]
    print("=== Выберите метод заполнения пропусков ===")
    for i, m in enumerate(possible_methods, start=1):
        print(f"{i}. {m}")
    choice = input("Введите номер (по умолчанию 1 = mean): ")
    try:
        idx = int(choice)
        if 1 <= idx <= len(possible_methods):
            return possible_methods[idx - 1]
    except ValueError:
        pass
    return "mean"


def detect_outliers(
    df: pd.DataFrame,
    column: str | None = None,
    *,
    method: Literal["iqr", "z-score", "mahalanobis", "isolation_forest"] = "iqr",
    threshold: float = 1.5,
) -> pd.DataFrame | pd.Series:
    """
    • `column is None` → DataFrame-маска (для IQR / Z-score)
    • `column` задан  → Series-маска для этого столбца

    Для Mahalanobis и Isolation-Forest всегда вычисляется по всем
    числовым столбцам, а результат при `column is None` сводится к Series,
    где *True* — строка-выброс.
    """
    if df.empty:
        empty = pd.DataFrame(False, index=df.index, columns=df.columns)
        return empty[column] if column else empty

    num_cols = df.select_dtypes(include=[np.number]).columns
    mask_df = pd.DataFrame(False, index=df.index, columns=num_cols)

    try:
        # --- IQR - Tukey -----------------------------------------------------------------
        if method == "iqr":
            cols = [column] if column else num_cols
            for col in cols:
                q1, q3 = df[col].quantile([0.25, 0.75])
                iqr = q3 - q1
                lb, ub = q1 - threshold * iqr, q3 + threshold * iqr
                mask_df[col] = (df[col] < lb) | (df[col] > ub)

        # --- Z-score ----------------------------------------------------------------------
        elif method == "z-score":
            cols = [column] if column else num_cols
            for col in cols:
                z = zscore(df[col].astype(float), nan_policy="omit")
                mask_df[col] = np.abs(z) > threshold

        # --- Mahalanobis ------------------------------------------------------------------
        elif method == "mahalanobis":
            num = df[num_cols].dropna()
            if len(num) >= 2:
                cov = EmpiricalCovariance().fit(num)
                dist = cov.mahalanobis(num)
                bad_idx = num.index[dist > threshold]
                mask_df.loc[bad_idx, :] = True

        # --- Isolation Forest -------------------------------------------------------------
        elif method == "isolation_forest":
            num = df[num_cols].dropna()
            if not num.empty:
                cont = threshold if 0 < threshold < 1 else 0.05
                iso = IsolationForest(contamination=cont, random_state=42)
                preds = iso.fit_predict(num)  # -1 = outlier
                bad_idx = num.index[preds == -1]
                mask_df.loc[bad_idx, :] = True

        else:
            raise ValueError(f"Unsupported method '{method}'")

    except Exception as exc:
        logger.exception("[detect_outliers] непредвидённая ошибка: %s", exc)

    # --- what to return ------------------------------------------------------------------
    if column:
        return mask_df[column]

    if method in {"iqr", "z-score"}:
        return mask_df

    # для Mahalanobis / IsolationForest вернём Series-по-строкам
    return mask_df.any(axis=1)


# --------------------------------------------------------------------------------------
#            — helpers —
# --------------------------------------------------------------------------------------
def _numeric_cols(df: pd.DataFrame, cols: list[str], allow_mixed: bool, ctx: str) -> list[str]:
    ok: list[str] = []
    for c in cols:
        if pd.api.types.is_numeric_dtype(df[c]):
            ok.append(c)
        elif not allow_mixed:
            logger.warning("[%s] '%s' не числовой – пропущен", ctx, c)
    return ok


def _fill_missing_with_ml(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    """
    Простейшая ML-импутация:
    • числовые → RandomForestRegressor
    • категориальные → LogisticRegression
    """
    out = df.copy()

    for col in cols:
        nan_mask = out[col].isna()
        if not nan_mask.any():
            continue

        X_other = out.loc[~nan_mask, out.columns.difference([col])]
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
            preds = model.predict(out.loc[nan_mask, X_other.columns])
            out.loc[nan_mask, col] = preds
        except Exception as exc:
            logger.warning("[_fill_missing_with_ml] %s – оставил NaN", exc)

    return out
