from __future__ import annotations

from typing import TYPE_CHECKING
import pandas as pd

if TYPE_CHECKING:
    # • Во время type-checking остаётся Series[bool] и т. д.
    BoolS = pd.Series[bool]
    NumS  = pd.Series[float | int]
else:
    # • Во время runtime — безопасный псевдоним
    BoolS = pd.Series
    NumS  = pd.Series
