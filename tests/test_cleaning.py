import pandas as pd
from genoscope.data_analysis.data_filtering import filter_outliers
from genoscope.data_analysis.data_cleaning import detect_outliers


def test_detect_outliers_mask():
    df = pd.DataFrame({"val": [1, 2, 3, 1000]})
    mask = detect_outliers(df, method="iqr")
    assert mask["val"].sum() == 1  # ровно один выброс


def test_filter_outliers_removed():
    df = pd.DataFrame({"val": [1, 2, 3, 1000]})
    df_clean = filter_outliers(df, "val", method="iqr")
    assert len(df_clean) == 3


def test_detect_outliers_iqr():
    df = pd.DataFrame({"x": [1, 2, 3, 1000]})
    mask = detect_outliers(df, column="x", method="iqr")
    assert mask.iloc[-1]  # 1000 — выброс
    assert mask.sum() == 1  # Только 1 выброс


def test_detect_outliers_zscore_large():
    values = [10] * 10 + [1000]
    df = pd.DataFrame({"x": values})
    mask = detect_outliers(df, column="x", method="z-score", threshold=2)
    assert mask.iloc[-1]
    assert mask.sum() == 1


def test_detect_outliers_mahalanobis_large():
    values = [10] * 10 + [1000]
    df = pd.DataFrame({"x": values, "y": values})
    mask = detect_outliers(df, method="mahalanobis", threshold=9.21)
    # mask — Series булевых значений по строкам
    assert mask.iloc[-1]
    assert mask.sum() == 1


def test_detect_outliers_isolation_forest():
    values = [10] * 15 + [1000]
    df = pd.DataFrame({"x": values, "y": values})
    mask = detect_outliers(df, method="isolation_forest", threshold=0.05)
    # Проверяем, что выброс выделился (из-за рандома их может быть больше одного)
    assert mask.iloc[-1]


def test_detect_outliers_no_outliers():
    df = pd.DataFrame({"x": [10, 12, 11, 13]})
    mask = detect_outliers(df, column="x", method="iqr")
    assert not mask.any()
