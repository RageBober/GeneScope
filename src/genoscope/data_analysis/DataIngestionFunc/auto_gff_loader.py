from __future__ import annotations

import os
from collections.abc import Generator
from collections.abc import Sequence

import pandas as pd

from .gff_enhanced import chunk_read_gff
from .gff_enhanced import load_gff_advanced


def auto_load_gff(  # ← полные аннотации
    file_path: str,
    *,
    size_threshold_mb: int = 200,
    chunk_size: int = 200_000,
    filter_types: Sequence[str] | None = None,
    parse_attrs: bool = True,
    dbfn: str = ":memory:",
    disk_mode: bool = True,
) -> pd.DataFrame | Generator[pd.DataFrame, None, None] | None:
    """
    Автоматически выбирает стратегию загрузки GFF/ GTF.

    • Если размер файла > ``size_threshold_mb`` ― стримим чанками;
    • иначе читаем целиком через :pyfunc:`load_gff_advanced`.

    Returns
    -------
    DataFrame | Generator[DataFrame, None, None]
    """
    size_mb: float = os.path.getsize(file_path) / (1024 * 1024)
    print(f"[auto_load_gff] Файл {file_path} ~ {size_mb:.2f} МБ")

    if size_mb > size_threshold_mb:
        print("[auto_load_gff] Файл большой → chunk_read_gff ...")
        return chunk_read_gff(
            file_path,
            chunk_size=chunk_size,
            filter_types=filter_types,
        )
    print("[auto_load_gff] Файл небольшой → load_gff_advanced ...")
    return load_gff_advanced(
        file_path,
        dbfn=dbfn,
        force=True,
        filter_types=filter_types,
        parse_attrs=parse_attrs,
        disk_mode=disk_mode,
    )
