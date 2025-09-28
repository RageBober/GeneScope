# tests/conftest.py
from __future__ import annotations

import json
import sys
from collections.abc import Iterator
from pathlib import Path

import pytest


# ────────────────────────────────────────────────────────────────────────
# существующие фикстуры
# ────────────────────────────────────────────────────────────────────────
@pytest.fixture(scope="session")
def data_dir() -> Path:
    """tests/assets, чтобы не писать длинные пути в каждом тесте."""
    return Path(__file__).parent / "assets"


@pytest.fixture
def csv_file(tmp_path: Path) -> Path:
    f = tmp_path / "sample.csv"
    f.write_text("a,b\n1,2\n3,4\n", encoding="utf-8")
    return f


@pytest.fixture
def json_file(tmp_path: Path) -> Path:
    f = tmp_path / "sample.json"
    json.dump({"id": 1, "val": 42}, f.open("w", encoding="utf-8"))
    return f


@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
    f = tmp_path / "sample.fasta"
    f.write_text(">seq1\nATGC\n>seq2\nGGCC\n", encoding="utf-8")
    return f


# ────────────────────────────────────────────────────────────────────────
# 👇 Новая фикстура-патч для cyvcf2.VCF (работает во всех тестах)
# ────────────────────────────────────────────────────────────────────────
class _DummyRecord:
    """Мини-эмуляция одной записи VCF из cyvcf2."""

    def __init__(self, chrom: str, pos: int) -> None:
        self.CHROM = chrom
        self.POS = pos
        self.ID = "."
        self.REF = "A"
        self.ALT = ["G"]
        self.QUAL = 60
        self.FILTER = "PASS"


class _DummyVCF:
    """
    Эмуляция класса cyvcf2.VCF: принимает path и итерируется
    по списку _DummyRecord.
    """

    def __init__(self, path: str) -> None:
        # Можно варьировать по path, если нужно несколько сценариев
        self._records = [_DummyRecord("1", 100), _DummyRecord("2", 200)]

    def __iter__(self) -> Iterator[_DummyRecord]:
        return iter(self._records)


@pytest.fixture(autouse=True)
def patch_cyvcf2(monkeypatch: pytest.MonkeyPatch) -> None:
    """
    Подменяет ``cyvcf2.VCF`` во всех тестах, чтобы исключить
    тяжёлую зависимость и реальные файлы.
    """
    fake_module = type("FakeCyvcf2", (), {"VCF": _DummyVCF})
    monkeypatch.setitem(sys.modules, "cyvcf2", fake_module)
