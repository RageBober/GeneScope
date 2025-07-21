# utils/logging_helpers.py
from __future__ import annotations

import logging
from functools import wraps
from typing import Callable, ParamSpec, TypeVar

P = ParamSpec("P")  # параметры оборачиваемой функции
R = TypeVar("R")  # возвращаемое значение

logger = logging.getLogger(__name__)


def log_exceptions(fn: Callable[P, R]) -> Callable[P, R]:
    """
    Декоратор. Если в обёрнутой функции возникает исключение —
    • пишет stack-trace через logger.exception,
    • повторно выбрасывает исключение, чтобы оно дошло до тестов/CI.
    """

    @wraps(fn)
    def _wrapper(*args: P.args, **kwargs: P.kwargs) -> R:  # type: ignore[override]
        try:
            return fn(*args, **kwargs)
        except Exception:  # noqa: BLE001 — ловим всё намеренно
            logger.exception("[%s] unhandled exception", fn.__name__)
            raise

    return _wrapper
