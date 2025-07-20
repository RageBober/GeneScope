import pandas as pd
from genoscope.data_analysis.data_ingestion import load_csv
import pytest

def test_load_csv_basic(data_dir):
    df = load_csv(data_dir / "mini.csv")
    assert isinstance(df, pd.DataFrame)
    assert len(df) == 3

# ---------- Краевые случаи -------------------------------------------------
@pytest.mark.parametrize("bad_path", ["no_such_file.csv", ""])
def test_load_csv_missing_file(tmp_path, bad_path):
    """Функция должна вызвать SystemExit, если файла нет или путь пустой."""
    with pytest.raises(SystemExit):
        load_csv(tmp_path / bad_path)          # <-- здесь используем фиктивный путь


def test_load_csv_empty(tmp_path):
    """Пустой файл → EmptyDataError внутри, функция вернёт None."""
    empty = tmp_path / "empty.csv"
    empty.write_text("")                       # создаём пустой файл

    df = load_csv(empty)
    assert df is None