import pandas as pd

from genoscope.data_analysis.data_filters.filter_ngs import filter_ngs


def test_filter_ngs_basic():
    df = pd.DataFrame({"READS": [5, 12, 30]})
    assert len(filter_ngs(df, min_reads=10)) == 2
