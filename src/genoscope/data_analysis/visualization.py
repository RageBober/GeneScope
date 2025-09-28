# data_analysis/visualization.py
from __future__ import annotations

import logging
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

logger = logging.getLogger(__name__)


def plot_correlation_matrix(
    df: pd.DataFrame, output_path: str | Path | None = None
) -> None:
    """
    Строит тепловую карту корреляционной матрицы DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
    output_path : str | Path | None
        Если указан — PNG сохраняется по этому пути.
    """
    try:
        corr = df.corr(numeric_only=True)
        plt.figure(figsize=(10, 8))
        sns.heatmap(corr, annot=True, cmap="coolwarm", fmt=".2f")
        plt.title("Correlation Matrix")
        if output_path:
            plt.savefig(Path(output_path).expanduser())
        if show_plot:

            plt.show()

        else:

            plt.close()
    except Exception as exc:
        logger.exception("plot_correlation_matrix: %s", exc)


def plot_distributions(df: pd.DataFrame, output_dir: str | Path | None = None) -> None:
    """
    Рисует гистограммы + KDE для всех числовых столбцов.

    Parameters
    ----------
    df : pd.DataFrame
    output_dir : str | Path | None
        Директория для PNG-файлов; если None — только вывод на экран.
    """
    try:
        num_cols = df.select_dtypes(include=["number"]).columns
        for col in num_cols:
            plt.figure()
            sns.histplot(df[col].dropna(), kde=True)
            plt.title(f"Distribution of {col}")
            if output_dir:
                out = Path(output_dir).expanduser() / f"{col}_distribution.png"
                plt.savefig(out)
            if show_plot:

                plt.show()

            else:

                plt.close()
    except Exception as exc:
        logger.exception("plot_distributions: %s", exc)


def plot_pca(
    pca_result: pd.DataFrame,
    target_column: Sequence[float | int] | pd.Series[Any] | None = None,
) -> None:
    """
    Визуализирует двумерный результат PCA.

    Parameters
    ----------
    pca_result : pd.DataFrame
        DataFrame с колонками «PC1», «PC2».
    target_column : Sequence | pd.Series | None
        Вектор цвета точек; если None — однотоннаяScatter-диаграмма.
    """
    try:
        plt.figure()
        plt.scatter(
            pca_result["PC1"],
            pca_result["PC2"],
            c=target_column if target_column is not None else "tab:blue",
            cmap="viridis",
        )
        plt.xlabel("PC1")
        plt.ylabel("PC2")
        plt.title("PCA Visualization")
        if target_column is not None:
            plt.colorbar(label="Target")
        if show_plot:

            plt.show()

        else:

            plt.close()
    except Exception as exc:
        logger.exception("plot_pca: %s", exc)
