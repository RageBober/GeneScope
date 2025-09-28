"""
Конфигурация для pytest тестирования.
"""

import shutil
import tempfile
from collections.abc import Generator
from pathlib import Path

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from genoscope.api.main import app


@pytest.fixture
def client() -> TestClient:
    """Тестовый клиент FastAPI."""
    return TestClient(app)


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """Временная директория для тестов."""
    temp_path = Path(tempfile.mkdtemp())
    try:
        yield temp_path
    finally:
        shutil.rmtree(temp_path, ignore_errors=True)


@pytest.fixture
def sample_csv_file(temp_dir: Path) -> Path:
    """Создает тестовый CSV файл."""
    csv_path = temp_dir / "test_data.csv"
    data = {
        "CHROM": ["1", "2", "3", "X"],
        "POS": [100, 200, 300, 400],
        "REF": ["A", "T", "G", "C"],
        "ALT": ["T", "C", "A", "G"],
        "GENE": ["GENE1", "GENE2", "GENE3", "GENE4"],
        "AF": [0.01, 0.05, 0.1, 0.5],
        "QUAL": [30, 40, 50, 60],
    }
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def sample_vcf_content() -> str:
    """Содержимое тестового VCF файла."""
    return """##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	T	30	PASS	AF=0.01
2	200	.	T	C	40	PASS	AF=0.05
3	300	.	G	A	50	PASS	AF=0.1
X	400	.	C	G	60	PASS	AF=0.5
"""


@pytest.fixture
def sample_dataframe() -> pd.DataFrame:
    """Тестовый DataFrame."""
    return pd.DataFrame(
        {
            "A": [1, 2, 3, 4, 5],
            "B": [10, 20, 30, 40, 50],
            "C": ["a", "b", "c", "d", "e"],
        }
    )


@pytest.fixture
def dataframe_with_missing() -> pd.DataFrame:
    """DataFrame с пропущенными значениями."""
    return pd.DataFrame(
        {
            "A": [1, None, 3, 4, None],
            "B": [10, 20, None, 40, 50],
            "C": ["a", None, "c", "d", "e"],
        }
    )


@pytest.fixture
def dataframe_with_duplicates() -> pd.DataFrame:
    """DataFrame с дубликатами."""
    return pd.DataFrame(
        {
            "sample_id": ["S1", "S1", "S2", "S3", "S3"],
            "score": [10, 15, 20, 25, 20],
            "value": ["a", "b", "c", "d", "e"],
        }
    )


@pytest.fixture
def dataframe_with_outliers() -> pd.DataFrame:
    """DataFrame с выбросами."""
    return pd.DataFrame(
        {
            "normal": [1, 2, 3, 4, 5],
            "with_outlier": [1, 2, 100, 4, 5],  # 100 - выброс
            "negative": [-1, -2, -3, -4, -5],
        }
    )
