#!/usr/bin/env python
"""
Простой тест для проверки работы тестовой системы
"""

def test_basic():
    """Базовый тест для проверки pytest"""
    assert 1 + 1 == 2

def test_import():
    """Тест импорта модулей"""
    try:
        import fastapi
        import sqlalchemy
        import pandas
        assert True
    except ImportError as e:
        assert False, f"Failed to import: {e}"
