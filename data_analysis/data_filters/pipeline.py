from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Callable, Dict, List
import logging
import pandas as pd

# низкоуровневые реализации лежат тут
from .. import data_filtering as _f 
"""data_filters.pipeline
~~~~~~~~~~~~~~~~~~~~~~~~
Высоко‑уровневый декларативный пайплайн фильтрации данных.

Позволяет собирать последовательность фильтров (`Step`) и применять
их к ``pandas.DataFrame`` одной строкой кода, не заботясь о внутренних
деталях импорта конкретных функций.

Библиотека строится поверх низко‑уровневых функций из ``data_analysis.data_filtering``.
Это избавляет от дублирования логики и сохраняет единое место правды.

Usage
-----
>>> from data_filters.pipeline import FilterPipeline
>>> pipe = (
...     FilterPipeline()
...     .add("percentile", column="depth", lower=1, upper=99)
...     .add("outliers", column="depth", method="iqr", threshold=1.5)
... )
>>> cleaned_df = pipe.run(df)

Либо конфиг‑ориентированный вариант::

    cfg = [
        {"type": "duplicates_filter", "subset": ["sample_id"], "keep": "first"},
        {"type": "nan_filter", "axis": 0, "how": "any"},
    ]
    cleaned = FilterPipeline.from_config(cfg).run(df)
"""




from data_analysis.data_filtering import (
    filter_by_multiple_conditions,
    filter_by_custom_function,
    filter_by_percentile,
    filter_outliers,
    filter_by_value_range,
    filter_by_nan,
    filter_duplicates,
    filter_by_correlation,
    filter_by_variance,
    filter_by_cluster,
)

logger: logging.Logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# 1. Реестр низко‑уровневых функций
# ---------------------------------------------------------------------------
_FilterFunc = Callable[[pd.DataFrame, Dict[str, Any]], pd.DataFrame]


def _wrap(func: Callable, *p_keys: str, **defaults: Any) -> _FilterFunc:
    """Преобразует функцию *func* в стандартный интерфейс ``(df, params) -> df``."""

    def _inner(df: pd.DataFrame, params: Dict[str, Any]) -> pd.DataFrame:  # type: ignore[name-defined]
        kwargs = {k: params.get(k, defaults.get(k)) for k in p_keys}
        return func(df, **kwargs)

    return _inner


_FILTER_REGISTRY: Dict[str, _FilterFunc] = {
    "multiple_conditions": _wrap(filter_by_multiple_conditions, "conditions", defaults={"conditions": []}),
    "custom_function": _wrap(filter_by_custom_function, "func"),
    "percentile": _wrap(filter_by_percentile, "column", "lower", "upper", defaults={"lower": 5, "upper": 95}),
    "outliers": _wrap(
        filter_outliers,
        "column",
        "method",
        "threshold",
        defaults={"method": "iqr", "threshold": 1.5},
    ),
    "value_range": _wrap(filter_by_value_range, "column", "min_value", "max_value"),
    "nan_filter": _wrap(
        filter_by_nan,
        "how",
        "axis",
        "subset",
        defaults={"how": "any", "axis": 0, "subset": None},
    ),
    "duplicates_filter": _wrap(filter_duplicates, "keep", "subset", defaults={"keep": "first"}),
    "correlation_filter": _wrap(filter_by_correlation, "corr_threshold", defaults={"corr_threshold": 0.95}),
    "variance_filter": _wrap(filter_by_variance, "var_threshold", defaults={"var_threshold": 0.0}),
    "cluster_filter": _wrap(
        filter_by_cluster,
        "method",
        "n_clusters",
        "drop_outliers",
        defaults={"method": "kmeans", "n_clusters": 2, "drop_outliers": True},
    ),
}

# ---------------------------------------------------------------------------
# 2. Dataclass Step — атомарная операция пайплайна
# ---------------------------------------------------------------------------


@dataclass
class Step:
    """Один шаг фильтрации."""

    type: str
    params: Dict[str, Any]

    def apply(self, df: pd.DataFrame) -> pd.DataFrame:
        if self.type not in _FILTER_REGISTRY:
            raise ValueError(f"[FilterPipeline] Неизвестный тип фильтра: '{self.type}'.")
        return _FILTER_REGISTRY[self.type](df, self.params)


# ---------------------------------------------------------------------------
# 3. Основной класс FilterPipeline
# ---------------------------------------------------------------------------


class FilterPipeline:
    """Гибкий последовательный пайплайн фильтрации ``pandas.DataFrame``."""

    def __init__(self, steps: List[Step] | None = None):
        self.steps: List[Step] = list(steps) if steps else []

    # -- construction ------------------------------------------------------
    def add(self, step_type: str, /, **params: Any) -> "FilterPipeline":
        """Добавить шаг и вернуть self (для цепочек)."""
        self.steps.append(Step(step_type, params))
        return self

    @classmethod
    def from_config(cls, cfg: List[Dict[str, Any]]) -> "FilterPipeline":
        """Создать пайплайн из списка словарей.

        Пример::
            cfg = [
                {"type": "duplicates_filter", "subset": ["id"]},
                {"type": "nan_filter", "axis": 0},
            ]
            pipe = FilterPipeline.from_config(cfg)
        """
        pipe = cls()
        for step_cfg in cfg:
            step_type = step_cfg.get("type")
            if not step_type:
                raise ValueError("Каждый шаг должен содержать ключ 'type'.")
            params = {k: v for k, v in step_cfg.items() if k != "type"}
            pipe.add(step_type, **params)
        return pipe

    # -- execution ---------------------------------------------------------
    def run(self, df: pd.DataFrame, *, inplace: bool = False, verbose: bool | int = False) -> pd.DataFrame:
        """Применить все шаги к DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Входные данные.
        inplace : bool, default False
            True → модифицируем ``df``; False → работаем на копии.
        verbose : bool | int, default False
            0 / False  → тихий режим,
            1 / True   → basic‑логи в stdout,
            2          → подробные логи (через ``logging``).
        """
        if not inplace:
            df = df.copy()

        for i, step in enumerate(self.steps, 1):
            if verbose:
                print(f"[FilterPipeline] ▶ step {i}/{len(self.steps)} — {step.type}")
            before = len(df)
            df = step.apply(df)
            after = len(df)
            if verbose == 2:
                logger.info("step %s → rows: %d → %d", step.type, before, after)
        return df

    # -- misc --------------------------------------------------------------
    def __repr__(self) -> str:  # pragma: no cover — косметика
        parts = ", ".join(step.type for step in self.steps) or "<empty>"
        return f"<FilterPipeline [{parts}]>"

    def __len__(self) -> int:
        return len(self.steps)

def run_filters(df: pd.DataFrame, steps: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Применяет к ``df`` последовательность фильтров, описанных в *steps*.

    Принимаются оба формата записи шага:

        {"name": "outliers", "params": {"column": "val"}}
        {"type": "outliers", "column": "val", ...}

    Первый («name» + «params») встречается в тестах-примерках.
    Второй — «плоский», совместим с ``apply_filter_pipeline`` напрямую.
    """
    normalized: list[dict[str, Any]] = []

    for s in steps:
        # поддержка обоих ключей: 'type' (новый) и 'name' (старый)
        step_type = s.get("type") or s.get("name")
        if not step_type:
            raise ValueError("Каждый шаг должен содержать ключ 'type' или 'name'.")

        # распаковываем вложенный словарь params (если есть)
        params: dict[str, Any] = s.get("params", {})
        normalized.append({"type": step_type, **params})

    # делегируем низкоуровневой функции-движку
    return _f.apply_filter_pipeline(df, normalized)

__all__ = [
    "FilterPipeline",
    "Step",
]
