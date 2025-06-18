# tests/test_variants_filter.py
import pandas as pd
from data_analysis.data_filters.filter_variants import filter_variants


def test_filter_variants_basic():
    df = pd.DataFrame({"QUAL": [10, 35, 60]})
    assert len(filter_variants(df, min_qual=30)) == 2
