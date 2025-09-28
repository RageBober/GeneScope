# main_large.py

import os
import pandas as pd
import numpy as np
from data_analysis.data_ingestion import load_data
from data_analysis.data_cleaning import remove_duplicates, handle_missing_values
from data_analysis.data_filtering import filter_by_multiple_conditions, filter_outliers
from data_analysis.analysis_core import extract_pca

print("Текущий рабочий каталог:", os.getcwd())


def generate_large_csv(file_path, target_size_mb=1024):
    """
    Генерирует CSV-файл примерно нужного размера.
    Здесь создается DataFrame с случайными значениями.

    target_size_mb - целевой размер файла в мегабайтах.
    """
    # Для приблизительной оценки: пусть каждый числовой элемент занимает ~8 байт,
    # плюс разделители и заголовки. Будем генерировать DataFrame с shape (n_rows, n_cols)
    # подбираем n_rows и n_cols так, чтобы итоговый файл был порядка target_size_mb.
    # Здесь выберем 120000 строк и 50 столбцов (примерно 120000*50*8 = 48 МБ чистых данных,
    # но с разделителями и заголовком получится больше – можно подстроить по необходимости).

    n_rows = 1200000
    n_cols = 500
    df = pd.DataFrame(
        np.random.rand(n_rows, n_cols), columns=[f"col{i}" for i in range(n_cols)]
    )
    df.to_csv(file_path, index=False)
    print(f"Large CSV file generated at {file_path} with shape {df.shape}.")


def main_large():
    file_path = "/home/asd/.virtualenvs/BioForge/1.GenoScope/data/large_sample.csv"
    # Если файла нет, создаем его
    if not os.path.exists(file_path):
        print("Большой файл не найден, генерируем...")
        generate_large_csv(file_path, target_size_mb=100)

    # Загружаем данные
    data = load_data(file_path, "csv")
    if data is None:
        print("Не удалось загрузить большой файл.")
        return
    print(f"Загружены данные: {data.shape}")

    # Применяем этапы обработки для проверки производительности:
    data = remove_duplicates(data)
    print("Дубликаты удалены.")

    data = handle_missing_values(data, method="mean")
    print("Пропущенные значения обработаны.")

    # Пример фильтрации: оставляем строки, где значение в первом столбце больше 0.5
    filtered_data = filter_by_multiple_conditions(data, ["col0 > 0.5"])
    print(f"После фильтрации по условиям: {filtered_data.shape}")

    # Пример фильтрации выбросов для одного столбца
    if "col0" in filtered_data.columns:
        filtered_data = filter_outliers(filtered_data, "col0", method="iqr")
        print(f"После фильтрации выбросов: {filtered_data.shape}")

    # Применяем PCA к числовым данным
    numeric_data = filtered_data.select_dtypes(include=["float64", "int"])
    try:
        pca_df = extract_pca(numeric_data, n_components=2)
        print(f"PCA результат: {pca_df.shape}")
    except Exception as e:
        print(f"Ошибка при выполнении PCA: {e}")


if __name__ == "__main__":
    main_large()
