"""
Конфигурация pytest для правильных путей
"""
import sys
from pathlib import Path

# Добавляем корневую директорию в PYTHONPATH
root_dir = Path(__file__).parent
sys.path.insert(0, str(root_dir))

# Конфигурация pytest
pytest_plugins = []
