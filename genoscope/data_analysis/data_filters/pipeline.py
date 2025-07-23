from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import (
    Any,
    Callable,
    Dict,
    List,
    Mapping,
    MutableMapping,
    Sequence,
    Literal,
    TypeAlias,   
)
import pandas as pd

from genoscope.data_analysis.data_filtering import (
    filter_by_cluster,
    filter_by_correlation,
    filter_by_custom_function,
    filter_by_multiple_conditions,
    filter_by_nan,
    filter_by_percentile,
    filter_by_value_range,
    filter_by_variance,
    filter_duplicates,
    filter_outliers,
)

# низкоуровневый движок
from .. import data_filtering as _f

# ──────────────────────────────────────────────────────────────────────────────
# Тип-алиасы
# ──────────────────────────────────────────────────────────────────────────────
logger = logging.getLogger(__name__)   
DF = pd.DataFrame
Params = Mapping[str, Any]
ParamsMutable = MutableMapping[str, Any]
FilterFunc = Callable[[DF, Params], DF]

# статический Literal: IDE молчит, mypy доволен
StepType: TypeAlias = Literal[
    "multiple_conditions",
    "custom_function",
    "percentile",
    "outliers",
    "value_range",
    "nan_filter",
    "duplicates_filter",
    "correlation_filter",
    "variance_filter",
    "cluster_filter",
]

# ──────────────────────────────────────────────────────────────────────────────
# Регистр низкоуровневых фильтров
# ──────────────────────────────────────────────────────────────────────────────
def _wrap(func: Callable[..., DF], *p_keys: str, defaults: Params | None = None) -> FilterFunc:
    defaults = defaults or {}

    def _inner(df: DF, params: Mapping[str, Any]) -> DF:
        kwargs = {k: params.get(k, defaults.get(k)) for k in p_keys}
        return func(df, **kwargs) 

    return _inner


_FILTER_REGISTRY: Dict[str, FilterFunc] = {
    "multiple_conditions": _wrap(
        filter_by_multiple_conditions,
        "conditions",
        defaults={"conditions": []},
    ),
    "custom_function": _wrap(filter_by_custom_function, "func"),
    "percentile": _wrap(
        filter_by_percentile,
        "column",
        "lower",
        "upper",
        defaults={"lower": 5, "upper": 95},
    ),
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
    "correlation_filter": _wrap(
        filter_by_correlation,
        "corr_threshold",
        defaults={"corr_threshold": 0.95},
    ),
    "variance_filter": _wrap(filter_by_variance, "var_threshold", defaults={"var_threshold": 0.0}),
    "cluster_filter": _wrap(
        filter_by_cluster,
        "method",
        "n_clusters",
        "drop_outliers",
        defaults={"method": "kmeans", "n_clusters": 2, "drop_outliers": True},
    ),
}

# ──────────────────────────────────────────────────────────────────────────────
# Шаг пайплайна
# ──────────────────────────────────────────────────────────────────────────────
@dataclass
class Step:
    name: StepType  # не «type» → нет A003
    params: Dict[str, Any]

    def apply(self, df: DF) -> DF:
        try:
            return _FILTER_REGISTRY[self.name](df, self.params)
        except KeyError as exc:
            raise ValueError(f"[FilterPipeline] Неизвестный фильтр «{self.name}».") from exc


# ──────────────────────────────────────────────────────────────────────────────
# Пайплайн
# ──────────────────────────────────────────────────────────────────────────────
class FilterPipeline:
    def __init__(self, steps: Sequence[Step] | None = None) -> None:
        self.steps: List[Step] = list(steps) if steps else []

    # — build —
    def add(self, step_name: StepType, /, **params: Any) -> "FilterPipeline":  # noqa: A002
        self.steps.append(Step(step_name, dict(params)))
        return self

    @classmethod
    def from_config(cls, cfg: Sequence[Dict[str, Any]]) -> "FilterPipeline":
        pipe = cls()
        for step in cfg:
            step_name = step.get("type") or step.get("name")
            if step_name is None:
                raise ValueError("Каждый шаг должен содержать ключ 'type'.")
            if step_name not in _FILTER_REGISTRY:
                raise ValueError(f"Неизвестный фильтр: «{step_name}».")
            params = {k: v for k, v in step.items() if k not in {"type", "name"}}
            pipe.add(step_name, **params) 
        return pipe

    # — run —
    def run(
        self,
        df: DF,
        *,
        inplace: bool = False,
        verbose: bool | int = False,
    ) -> DF:
        if not inplace:
            df = df.copy()

        for idx, step in enumerate(self.steps, 1):
            if verbose:
                print(f"[FilterPipeline] ▶ step {idx}/{len(self.steps)} — {step.name}")
            before = len(df)
            df = step.apply(df)
            after = len(df)
            if verbose == 2:
                logger.info(  # noqa: PLE1202
                    "step %s → rows: %d → %d",
                    step.name,
                    before,
                    after,
                )
        return df

    # misc
    def __repr__(self) -> str:
        return f"<FilterPipeline [{', '.join(s.name for s in self.steps) or '<empty>'}]>"

    def __len__(self) -> int:
        return len(self.steps)


# ──────────────────────────────────────────────────────────────────────────────
# Обёртка run_filters
# ──────────────────────────────────────────────────────────────────────────────
def run_filters(df: DF, steps: Sequence[Dict[str, Any]]) -> DF:
    normalized: list[dict[str, Any]] = []
    for raw in steps:
        step_name = raw.get("type") or raw.get("name")
        if step_name is None:
            raise ValueError("Каждый шаг должен иметь ключ 'type' или 'name'.")
        params = dict(raw.get("params", {}))
        normalized.append({"type": step_name, **params})
    return _f.apply_filter_pipeline(df, normalized) 


__all__: list[str] = ["FilterPipeline", "Step"]
