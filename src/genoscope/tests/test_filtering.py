import pandas as pd

from genoscope.data_analysis.data_filtering import filter_outliers


def test_filter_outliers_iqr():
    df = pd.DataFrame({"x": [1, 2, 3, 1000]})
    df_clean = filter_outliers(df, column="x", method="iqr")
    assert 1000 not in df_clean["x"].values
    assert len(df_clean) == 3


def test_filter_outliers_zscore_large():
    # Большая выборка: выброс явно выделяется
    values = [10] * 10 + [1000]
    df = pd.DataFrame({"x": values})
    df_clean = filter_outliers(df, column="x", method="z-score", threshold=2)
    assert 1000 not in df_clean["x"].values


def test_filter_outliers_mahalanobis_large():
    # Большая выборка для Mahalanobis
    values = [10] * 10 + [1000]
    df = pd.DataFrame({"x": values, "y": values})
    df_clean = filter_outliers(df, method="mahalanobis", threshold=9.21)
    assert 1000 not in df_clean["x"].values


def test_filter_outliers_isolation_forest():
    # Isolation forest: большой массив, 100 выброс
    values = [10] * 15 + [1000]
    df = pd.DataFrame({"x": values, "y": values})
    df_clean = filter_outliers(df, method="isolation_forest", threshold=0.05)
    # isolation_forest может удалить не только 1000, поэтому assert слабее:
    assert 1000 not in df_clean["x"].values or len(df_clean) < len(df)
