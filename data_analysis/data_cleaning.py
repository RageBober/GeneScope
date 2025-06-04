# data_cleaning.py
"""
Модуль для очистки и предварительной обработки данных GenoScope.
Включает функции удаления дубликатов, обработки пропусков, обнаружения выбросов и вспомогательные утилиты.
"""

import pandas as pd
import logging

from utils.outlier_utils import get_outlier_mask
from MLmodel.data_cleaning_AI.ml_imputation import fill_missing_with_ml

logger = logging.getLogger(__name__)


def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Удаляет дубликаты строк из DataFrame.

    Параметры:
        df (pd.DataFrame): Исходный DataFrame.

    Возвращает:
        pd.DataFrame: DataFrame без дубликатов.
    """
    if df.empty:
        logger.warning("[remove_duplicates] Пустой DataFrame → без изменений.")
        return df

    before_count = len(df)
    df_clean = df.drop_duplicates()
    after_count = len(df_clean)
    removed = before_count - after_count
    if removed > 0:
        logger.info(f"[remove_duplicates] Удалено дубликатов: {removed}")
    return df_clean


def handle_missing_values(
    df: pd.DataFrame,
    method: str = "ffill",
    columns: list = None,
    *,
    allow_mixed_types: bool = False,
) -> pd.DataFrame:
    """
    Обрабатывает пропущенные значения в DataFrame выбранным методом.

    Поддерживаются методы:
    'ffill', 'bfill', 'mean', 'median', 'mode', 'interpolate', 'ml', 'dialog'.

    Параметры:
        df (pd.DataFrame): Исходный DataFrame.
        method (str): Метод заполнения ('ffill', 'mean', 'ml' и т.д.).
        columns (list, optional): Список столбцов, где заполнять пропуски (по умолчанию все).
        allow_mixed_types (bool): Разрешать ли нечисловые столбцы для mean/median/interpolate.

    Возвращает:
        pd.DataFrame: Новый DataFrame с заполненными пропусками.
    """
    if df.empty:
        logger.warning(
            "[handle_missing_values] Пустой DataFrame, возврат без изменений."
        )
        return df

    if method == "dialog":
        # Пример псевдо-GUI, где пользователь выбирает метод в консоли
        chosen = ask_user_method()  # Эта функция ниже — просто пример
        logger.info(
            f"[handle_missing_values] Пользователь выбрал метод через dialog: {chosen}"
        )
        method = chosen

    if columns is None:
        columns = df.columns.tolist()

    # Проверим, что все указанные столбцы существуют
    missing_cols = [c for c in columns if c not in df.columns]
    if missing_cols:
        logger.warning(
            f"[handle_missing_values] Столбцы {missing_cols} не найдены. Пропускаем их."
        )
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
            num_cols = _check_numeric_columns(
                df_copy, columns, allow_mixed_types, context="mean"
            )
            mean_vals = df_copy[num_cols].mean(numeric_only=True)
            df_copy[num_cols] = df_copy[num_cols].fillna(mean_vals)

        elif method == "median":
            # Медиана — только числовые столбцы
            num_cols = _check_numeric_columns(
                df_copy, columns, allow_mixed_types, context="median"
            )
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
    df: pd.DataFrame, columns: list[str], allow_mixed: bool, context: str
) -> list[str]:
    """
    Проверяет список столбцов на числовой тип.

    Параметры:
        df (pd.DataFrame): Данные.
        columns (list): Список столбцов.
        allow_mixed (bool): Разрешать ли нечисловые столбцы.
        context (str): Контекст (для логирования).

    Возвращает:
        list: Список числовых столбцов.
    """
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
    Имитирует диалоговый выбор метода заполнения пропусков (например, через консоль).

    Возвращает:
        str: Выбранный пользователем метод заполнения.
    """
    # В реальности здесь был бы вызов:
    # - Окно Tkinter
    # - radio buttons или dropdown
    # - кнопка OK
    # Но для примера просто спросим в консоли:
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
    method: str = "iqr",
    threshold: float = 1.5,
) -> pd.Series | pd.DataFrame:
    """
    Обнаружение выбросов в DataFrame.

    Возвращает булеву маску, где выбросы обозначены как True.

    Параметры:
        df (pd.DataFrame): Данные для анализа.
        column (str или None): Имя столбца или None (по всем числовым).
        method (str): Метод ('iqr', 'z-score', 'mahalanobis', 'isolation_forest').
        threshold (float): Пороговый параметр.

    Возвращает:
        pd.Series или pd.DataFrame: Маска выбросов.
    """
    return get_outlier_mask(df, column, method, threshold)
