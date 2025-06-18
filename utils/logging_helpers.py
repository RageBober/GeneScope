# utils/logging_helpers.py
from functools import wraps
import logging

logger = logging.getLogger(__name__)

def log_exceptions(fn):
    """
    Оборачивает функцию; если внутри возникает исключение —
    ▸ пишет stack trace в лог через logger.exception  
    ▸ пробрасывает исключение дальше (чтобы тесты/CI не «замалчивались»)
    """
    @wraps(fn)
    def _wrapper(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except Exception:          # ловим ВСЁ, логируем 1 строкой + traceback
            logger.exception("[%s] unhandled exception", fn.__name__)
            raise
    return _wrapper
