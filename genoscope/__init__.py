"""Top-level package for GenoScope."""

from . import data_analysis, mlmodel, utils

__all__ = ["data_analysis", "mlmodel", "utils"]

"""
Обратная совместимость: временно даёт писать

    import data_analysis.<submod>

а на деле перенаправляет вызов в новый пакет

    genoscope.data_analysis.<submod>

Как только весь код будет переведён на новую схему импорта,
этот файл можно удалить.
"""
from importlib import import_module as _imp
import sys as _sys

# Подключаем настоящий пакет
_real_pkg = _imp("genoscope.data_analysis")

# Регистрируем его под старым именем, чтобы `import data_analysis...` работал
_sys.modules[__name__] = _real_pkg

# Экспортируем все атрибуты (dir(), help(), mypy и т. д.)
globals().update(_real_pkg.__dict__)

# Чистим namespace вспомогательных переменных
del _imp, _sys, _real_pkg
