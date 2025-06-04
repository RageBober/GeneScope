from data_analysis.data_ingestion import load_data
from data_analysis.data_cleaning import remove_duplicates, handle_missing_values
from data_analysis.data_filtering import (
    filter_by_multiple_conditions,
    filter_by_custom_function,
    filter_by_percentile,
    filter_outliers,
)
from data_analysis.analysis_core import extract_pca
from data_analysis.visualization import plot_pca
from sklearn.preprocessing import LabelEncoder
import pandas as pd
import tkinter as tk
from interface import GenoScopeApp
import os

print(os.getcwd())


def apply_filters(data: pd.DataFrame) -> pd.DataFrame:
    """
    Применение различных фильтров к данным.

    Параметры:
        data (pd.DataFrame): Исходный DataFrame.

    Возвращает:
        pd.DataFrame: Отфильтрованные данные.
    """
    print("Применение фильтрации данных...")
    try:
        # Фильтрация по нескольким условиям
        filtered_data = filter_by_multiple_conditions(data, ["A > 10", "B < 50"])
        print(f"Данные после фильтрации по условиям:\n{filtered_data}")

        # Фильтрация с пользовательской функцией
        filtered_data = filter_by_custom_function(
            filtered_data, lambda row: row["A"] % 2 == 0 and row["B"] > 20
        )
        print(f"Данные после пользовательской фильтрации:\n{filtered_data}")

        # Фильтрация по персентилям
        if "A" in filtered_data.columns:
            filtered_data = filter_by_percentile(
                filtered_data, "A", lower_percentile=10, upper_percentile=90
            )
            print(f"Данные после фильтрации по персентилям:\n{filtered_data}")

        # Фильтрация выбросов
        if "A" in filtered_data.columns:
            filtered_data = filter_outliers(filtered_data, "A", method="iqr")
            print(f"Данные после фильтрации выбросов:\n{filtered_data}")

        return filtered_data
    except Exception as e:
        print(f"Ошибка при фильтрации данных: {e}")
        return data


def perform_pca(data: pd.DataFrame) -> pd.DataFrame:
    """
    Выполнение анализа PCA на числовых данных.

    Параметры:
        data (pd.DataFrame): Исходный DataFrame.

    Возвращает:
        pd.DataFrame: DataFrame с главными компонентами.
    """
    # Если имеется целевой столбец, исключаем его из признаков
    if "target_column" in data.columns:
        target_column = data["target_column"]
        df_features = data.drop(columns=["target_column"])
    else:
        target_column = None
        df_features = data

    try:
        # Выполнение PCA с использованием функции из analysis_core
        pca_df = extract_pca(df_features, n_components=2)
        print(f"PCA Result:\n{pca_df}")

        # Визуализация результатов PCA
        plot_pca(pca_df, target_column)
        return pca_df
    except Exception as e:
        print(f"Ошибка PCA анализа: {e}")
        return pd.DataFrame()


def main():
    """
    Основной pipeline для GenoScope.
    """
    # Шаг 1: Загрузка данных
    file_path = "data/sample.csv"  # Укажите путь к вашему файлу
    file_type = "csv"  # Укажите тип файла
    data = load_data(file_path, file_type)

    if data is None:
        print("Не удалось загрузить данные.")
        return

    print(f"Загруженные данные:\n{data.head()}")

    # Шаг 2: Очистка данных
    data = remove_duplicates(data)
    print("Удалены дубликаты.")

    data = handle_missing_values(data, method="mean")
    print("Обработаны пропущенные значения.")

    # Кодирование категориальных данных
    for column in data.select_dtypes(include=["object", "category"]).columns:
        le = LabelEncoder()
        data[column] = le.fit_transform(data[column])
    print("Кодированы категориальные данные (если были).")

    # Шаг 3: Фильтрация данных
    data = apply_filters(data)

    # Шаг 4: PCA
    numeric_data = data.select_dtypes(include=["float64", "int"])
    pca_result = perform_pca(numeric_data)

    # Шаг 5: Вывод результатов
    print(f"PCA Result:\n{pca_result}")


if __name__ == "__main__":
    # Запуск графического интерфейса
    root = tk.Tk()
    app = GenoScopeApp(root)
    root.mainloop()
