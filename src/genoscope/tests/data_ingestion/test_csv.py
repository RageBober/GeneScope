# tests/data_ingestion/test_csv.py
from __future__ import annotations

from collections.abc import Iterator
from pathlib import Path
from typing import cast

import pandas as pd
import pytest

from genoscope.data_analysis.data_ingestion import load_csv


@pytest.fixture
def tmp_csv(tmp_path: Path) -> Path:
    """Создаём временный CSV с простыми числовыми данными."""
    path = tmp_path / "sample.csv"
    df = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    df.to_csv(path, index=False)
    return path


def test_load_csv_plain(tmp_csv: Path) -> None:
    df = load_csv(str(tmp_csv))
    assert isinstance(df, pd.DataFrame)
    pd.testing.assert_frame_equal(  # type: ignore[arg-type]
        df.reset_index(drop=True),
        pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]}),
    )


def test_load_csv_chunked(tmp_csv: Path) -> None:
    chunks = load_csv(str(tmp_csv), chunksize=2)
    chunks = cast("Iterator[pd.DataFrame]", chunks)

    collected = pd.concat(list(chunks)).reset_index(drop=True)
    expected = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    pd.testing.assert_frame_equal(collected, expected)


def test_load_csv_empty(tmp_path: Path) -> None:
    empty_path = tmp_path / "empty.csv"
    empty_path.write_text("")  # создаём пустой файл

    assert load_csv(str(empty_path)) is None
