"""data_analysis.utils.columns
---------------------------
Утилиты для работы с колонками DataFrame.
"""

from __future__ import annotations

import logging
from typing import Sequence

import pandas as pd

logger = logging.getLogger(__name__)

def numeric(df: pd.DataFrame, cols: Sequence[str], allow_mixed: bool, ctx: str) -> list[str]:
    numeric_cols = []
    for col in cols:
        if col not in df.columns:
            continue
        if pd.api.types.is_numeric_dtype(df[col]):
            numeric_cols.append(col)
        else:
            if allow_mixed:
                logger.info(f"[numeric] '{col}' нечисловой, но allow_mixed_types=True => пропускаем предупреждение.")
            else:
                logger.warning(f"[numeric] '{col}' нечисловой, исключаем для метода '{ctx}'.")
    return numeric_cols