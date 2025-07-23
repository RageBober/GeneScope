import pandas as pd

from genoscope.data_analysis.data_ingestion import load_csv


def make_csv(path, rows):
    pd.DataFrame(rows).to_csv(path, index=False)


def test_load_csv_ok(tmp_path):
    csv_file = tmp_path / "good.csv"
    make_csv(csv_file, [{"a": 1, "b": 2}, {"a": 3, "b": 4}])

    df = load_csv(str(csv_file))
    assert df is not None
    assert len(df) == 2
    assert list(df.columns) == ["a", "b"]

    # optional: chunked reading
    chunks = list(load_csv(str(csv_file), chunksize=1))
    assert all(len(chunk) == 1 for chunk in chunks)


def test_load_csv_empty(tmp_path, caplog):
    empty_csv = tmp_path / "empty.csv"
    empty_csv.write_text("")

    with caplog.at_level("WARNING"):
        df = load_csv(str(empty_csv))
        assert df is None
        assert any("contains no data" in rec.message or "Loaded data is empty" in rec.message for rec in caplog.records)
