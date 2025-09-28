"""Базовые тесты GenoScope."""
import sys
from pathlib import Path

# Добавляем src в путь
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def test_imports():
    """Тест основных импортов."""
    from genoscope.main import GenoScopeProcessor
    from genoscope.core.validation import DataValidator
    # If we get here without ImportError, the test passes
    assert True

def test_processor_creation():
    """Тест создания процессора."""
    from genoscope.main import GenoScopeProcessor
    processor = GenoScopeProcessor()
    assert processor is not None

def test_validation():
    """Тест валидации."""
    from genoscope.core.validation import DataValidator
    valid, msg = DataValidator.validate_file_path("nonexistent.csv")
    assert isinstance(valid, bool)
    assert isinstance(msg, str)

if __name__ == "__main__":
    print("🧪 Запуск базовых тестов...")
    
    tests = [test_imports, test_processor_creation, test_validation]
    passed = 0
    
    for test in tests:
        if test():
            passed += 1
    
    print(f"\n📊 Результат: {passed}/{len(tests)} тестов прошло")
