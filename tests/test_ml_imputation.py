import pandas as pd
from MLmodel.data_cleaning_AI.ml_imputation import fill_missing_with_ml


def test_fill_missing_numeric():
    df = pd.DataFrame({"x": [1, 2, None, 4], "y": [5, 6, 7, 8]})
    result = fill_missing_with_ml(df, columns=["x"])
    assert not result["x"].isnull().any()


def test_fill_missing_categorical():
    df = pd.DataFrame({"cat": ["a", "b", None, "b"], "num": [1, 2, 3, 4]})
    result = fill_missing_with_ml(df, columns=["cat"])
    assert not result["cat"].isnull().any()
