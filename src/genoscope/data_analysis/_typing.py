from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

import pandas as pd

if TYPE_CHECKING:
    # Type-safe aliases for type checking
    BoolS = pd.Series[bool]
    NumS = pd.Series[float | int]
    DF: TypeAlias = pd.DataFrame[any]  # Proper type alias
    BoolSeries: TypeAlias = pd.Series[bool] 
    BoolFrame: TypeAlias = pd.DataFrame[bool]
else:
    # Runtime-safe aliases
    BoolS = pd.Series
    NumS = pd.Series
    DF = pd.DataFrame
    BoolSeries = pd.Series
    BoolFrame = pd.DataFrame

# Export the aliases
__all__ = ["BoolS", "NumS", "DF", "BoolSeries", "BoolFrame"]
