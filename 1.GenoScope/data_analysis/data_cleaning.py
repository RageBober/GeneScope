# data_cleaning.py
import pandas as pd
import numpy as np
import logging

from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.exceptions import NotFittedError
from scipy.stats import zscore

from MLmodel.data_cleaning_AI.ml_imputation import fill_missing_with_ml

logger = logging.getLogger(__name__)

def remove_duplicates(df: pd.DataFrame) -> pd.DataFrame:
    """
    Удаляет дубликаты из DataFrame.
    """
    if df.empty:
        logger.warning("[remove_duplicates] Поступил пустой DataFrame, возвращаем без изменений.")
        return df

    before_count = len(df)
    df_clean = df.drop_duplicates()
    after_count = len(df_clean)
    removed = before_count - after_count
    if removed > 0:
        logger.info(f"[remove_duplicates] Удалено дубликатов: {removed}")
    return df_clean


def handle_missing_values(df: pd.DataFrame, method: str = 'ffill',
    columns: list = None, *, allow_mixed_types: bool = False
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

    if method == 'dialog':
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
            logger.warning("[handle_missing_values] Ни одного подходящего столбца не осталось, возврат DF без изменений.")
            return df

    df_copy = df.copy()

    try:
        if method in ['ffill', 'bfill']:
            # Используем встроенные методы pandas fillna
            df_copy[columns] = df_copy[columns].fillna(method=method)

        elif method == 'mean':
            # Для среднего берем только числовые столбцы
            num_cols = _check_numeric_columns(df_copy, columns, allow_mixed_types, context='mean')
            mean_vals = df_copy[num_cols].mean(numeric_only=True)
            df_copy[num_cols] = df_copy[num_cols].fillna(mean_vals)

        elif method == 'median':
            # Медиана — только числовые столбцы
            num_cols = _check_numeric_columns(df_copy, columns, allow_mixed_types, context='median')
            median_vals = df_copy[num_cols].median(numeric_only=True)
            df_copy[num_cols] = df_copy[num_cols].fillna(median_vals)

        elif method == 'mode':
            # Моду вычисляем для каждого столбца отдельно
            # (для числовых и категориальных без ограничений).
            for col in columns:
                series_mode = df_copy[col].mode(dropna=True)
                if not series_mode.empty:
                    # Если мод несколько, берем первый
                    fill_val = series_mode.iloc[0]
                    df_copy[col] = df_copy[col].fillna(fill_val)

        elif method == 'interpolate':
            # Линейная интерполяция только для числовых
            num_cols = _check_numeric_columns(df_copy, columns, allow_mixed_types, context='interpolate')
            df_copy[num_cols] = df_copy[num_cols].interpolate(method='linear')

        elif method == 'ml':
            df_copy = fill_missing_with_ml(df_copy, columns)

        else:
            raise ValueError(f"[handle_missing_values] Unsupported method: '{method}'")

    except Exception as exc:
        logger.error(f"[handle_missing_values] Ошибка при заполнении методом '{method}': {exc}", exc_info=True)
        return df  # Возвращаем исходный DataFrame при ошибке

    # Логируем, сколько пропусков осталось
    leftover_nans = df_copy[columns].isna().sum().sum()
    logger.info(f"[handle_missing_values] Метод '{method}' завершён. Осталось пропусков: {leftover_nans}")
    return df_copy


# -------------------------------------------------------------------------
# Вспомогательная функция для проверки "числовых" столбцов
# -------------------------------------------------------------------------
def _check_numeric_columns(df: pd.DataFrame, columns: list, allow_mixed: bool, context: str) -> list:
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
                logger.warning(f"[_check_numeric_columns] '{col}' нечисловой, исключаем для метода '{context}'.")
    return numeric_cols


# -------------------------------------------------------------------------
# Пример ПСЕВДО-функции, как будто открываем GUI
# -------------------------------------------------------------------------
def ask_user_method() -> str:
    """
    Имитирует GUI-спрос у пользователя (в реальном коде это мог бы быть диалог Tkinter).
    Возвращает строку: 'mean' / 'median' / 'mode' / 'ml' / ...
    """
    # В реальности здесь был бы вызов:
    # - Окно Tkinter
    # - radio buttons или dropdown
    # - кнопка OK
    # Но для примера просто спросим в консоли:
    possible_methods = ['mean', 'median', 'mode', 'ml', 'interpolate', 'ffill', 'bfill']
    print("=== Выберите метод заполнения пропусков ===")
    for i, m in enumerate(possible_methods, start=1):
        print(f"{i}. {m}")
    choice = input("Введите номер (по умолчанию 1 = mean): ")
    try:
        idx = int(choice)
        if 1 <= idx <= len(possible_methods):
            return possible_methods[idx-1]
    except ValueError:
        pass
    return 'mean'

def detect_outliers(df: pd.DataFrame, method: str = 'iqr') -> pd.DataFrame:
    """
    Возвращает булев DataFrame того же размера: True=выброс.
    """
    if df.empty:
        logger.warning("[detect_outliers] Пустой DataFrame, возвращаем пустую mask-таблицу.")
        return pd.DataFrame(False, index=df.index, columns=df.columns)

    numeric_cols = df.select_dtypes(include=[np.number]).columns
    mask_outliers = pd.DataFrame(False, index=df.index, columns=df.columns)

    try:
        if method == 'iqr':
            for col in numeric_cols:
                q1 = df[col].quantile(0.25)
                q3 = df[col].quantile(0.75)
                iqr = q3 - q1
                lb = q1 - 1.5 * iqr
                ub = q3 + 1.5 * iqr
                mask_outliers[col] = (df[col] < lb) | (df[col] > ub)

        elif method == 'z-score':
            z_scores = df[numeric_cols].apply(zscore)
            outliers_bool = (z_scores.abs() > 3)
            for col in numeric_cols:
                mask_outliers[col] = outliers_bool[col]
        else:
            raise ValueError(f"[detect_outliers] Unsupported method: {method}")
    except Exception as e:
        logger.error(f"[detect_outliers] Ошибка: {e}", exc_info=True)

    count_outliers = mask_outliers.values.sum()
    logger.info(f"[detect_outliers] Метод '{method}', всего выбросов: {count_outliers}.")
    return mask_outliers