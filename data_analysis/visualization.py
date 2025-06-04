import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_correlation_matrix(df: pd.DataFrame, output_path: str | None = None) -> None:
    """
    Строит и отображает корреляционную матрицу для всех числовых признаков DataFrame.

    Параметры:
        df (pd.DataFrame): Входные данные.
        output_path (str, optional): Путь для сохранения изображения. Если None, только показывается на экране.

    Возвращает:
        None
    """
    try:
        corr = df.corr()
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap="coolwarm")
        plt.title("Correlation Matrix")
        if output_path:
            plt.savefig(output_path)
        plt.show()
    except Exception as e:
        print(f"Error plotting correlation matrix: {e}")


def plot_distributions(df: pd.DataFrame, output_dir: str | None = None) -> None:
    """
    Строит и отображает гистограммы распределения для всех числовых признаков DataFrame.

    Параметры:
        df (pd.DataFrame): Входные данные.
        output_dir (str, optional): Каталог для сохранения картинок. Если None, только показываются на экране.

    Возвращает:
        None
    """
    try:
        numeric_cols = df.select_dtypes(include=["float64", "int"]).columns
        for col in numeric_cols:
            plt.figure()
            sns.histplot(df[col], kde=True)
            plt.title(f"Distribution of {col}")
            if output_dir:
                plt.savefig(f"{output_dir}/{col}_distribution.png")
            plt.show()
    except Exception as e:
        print(f"Error plotting distributions: {e}")


def plot_pca(pca_result: pd.DataFrame, target_column: list | None = None) -> None:
    """
    Визуализирует результаты анализа главных компонент (PCA).

    Параметры:
        pca_result (pd.DataFrame): DataFrame с колонками 'PC1' и 'PC2' (и, опционально, целевой переменной).
        target_column (array-like, optional): Классы/метки для раскраски точек.

    Возвращает:
        None
    """
    try:
        plt.scatter(
            pca_result["PC1"], pca_result["PC2"], c=target_column, cmap="viridis"
        )
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("PCA Visualization")
        if target_column is not None:
            plt.colorbar(label="Target")
        plt.show()
    except Exception as e:
        print(f"Error plotting PCA: {e}")
