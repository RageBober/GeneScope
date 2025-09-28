import numpy as np
import pandas as pd

from genoscope.data_analysis.data_filters.pipeline import run_filters


def test_pipeline_outliers():
    df = pd.DataFrame({"val": np.r_[np.random.randn(10), 999]})
    steps = [{"name": "outliers", "params": {"column": "val"}}]
    df_clean = run_filters(df, steps)
    assert df_clean["val"].max() < 500  # выброс удалён
