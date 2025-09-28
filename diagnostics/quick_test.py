#!/usr/bin/env python3
"""
Quick Test Script for BioForge
Быстрая проверка основных проблем
"""

import sys
from pathlib import Path

# Добавляем src в путь для импортов
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))

def test_imports():
    """Тест импортов основных модулей"""
    print("\n🔗 Тестирование импортов...")
    
    results = []
    
    # Тест 1: main.py
    try:
        from genoscope.main import GenoScopeProcessor
        results.append("✅ main.py импортируется")
    except ImportError as e:
        results.append(f"❌ main.py: {e}")
    
    # Тест 2: interface.py
    try:
        from genoscope.interface import GenoScopeApp
        results.append("✅ interface.py импортируется")
    except ImportError as e:
        results.append(f"❌ interface.py: {e}")
    
    # Тест 3: data_ingestion.py
    try:
        from genoscope.data_analysis.data_ingestion import load_data
        results.append("✅ data_ingestion.py импортируется")
    except ImportError as e:
        results.append(f"❌ data_ingestion.py: {e}")
    
    # Тест 4: analysis_core.py
    try:
        from genoscope.data_analysis.analysis_core import extract_pca
        results.append("✅ analysis_core.py импортируется")
    except ImportError as e:
        results.append(f"❌ analysis_core.py: {e}")
    
    # Тест 5: visualization.py
    try:
        from genoscope.data_analysis.visualization import plot_correlation_matrix
        results.append("✅ visualization.py импортируется")
    except ImportError as e:
        results.append(f"❌ visualization.py: {e}")
    
    # Тест 6: API
    try:
        from genoscope.api.main import app
        results.append("✅ API модуль импортируется")
    except ImportError as e:
        results.append(f"❌ API модуль: {e}")
    
    # Тест 7: Pipeline
    try:
        from genoscope.pipeline.qc import QCMetrics
        from genoscope.pipeline.variant_calling import VariantStats
        results.append("✅ Pipeline модули импортируются")
    except ImportError as e:
        results.append(f"❌ Pipeline модули: {e}")
    
    return results

def test_file_structure():
    """Проверка структуры файлов"""
    print("\n📁 Проверка структуры проекта...")
    
    results = []
    
    required_files = [
        "src/genoscope/__init__.py",
        "src/genoscope/main.py",
        "src/genoscope/core/logging_config.py",
        "src/genoscope/data_analysis/__init__.py",
        "src/genoscope/data_analysis/data_ingestion.py",
        "src/genoscope/api/main.py",
        "requirements.txt",
        "pyproject.toml"
    ]
    
    for file_path in required_files:
        full_path = project_root / file_path
        if full_path.exists():
            results.append(f"✅ {file_path}")
        else:
            results.append(f"❌ {file_path} - НЕ НАЙДЕН")
    
    return results

def test_critical_issues():
    """Проверка критических проблем"""
    print("\n⚠️ Проверка критических проблем...")
    
    results = []
    
    # Проверка interface.py
    interface_path = project_root / "src" / "genoscope" / "interface.py"
    if not interface_path.exists():
        results.append("❌ КРИТИЧНО: Отсутствует interface.py (GUI модуль)")
    else:
        results.append("✅ interface.py существует")
    
    # Проверка валидации
    validation_path = project_root / "src" / "genoscope" / "core" / "validation.py"
    if not validation_path.exists():
        results.append("❌ КРИТИЧНО: Отсутствует модуль валидации файлов")
    else:
        results.append("✅ Модуль валидации существует")
    
    # Проверка логирования
    logging_path = project_root / "src" / "genoscope" / "core" / "logging_config.py"
    if logging_path.exists():
        content = logging_path.read_text()
        if '"handlers": ["console", "console_error"]' in content:
            results.append("⚠️ Проблема: Дублирование в логировании")
        else:
            results.append("✅ Логирование настроено корректно")
    
    # Проверка типов
    ingestion_path = project_root / "src" / "genoscope" / "data_analysis" / "data_ingestion.py"
    if ingestion_path.exists():
        content = ingestion_path.read_text()
        if '| None' in content:
            results.append("⚠️ Проблема: Использование Python 3.10+ синтаксиса типов")
        else:
            results.append("✅ Типы совместимы")
    
    return results

def main():
    print("=" * 60)
    print("🧬 BioForge Quick Test")
    print("=" * 60)
    
    # Тест структуры
    structure_results = test_file_structure()
    for result in structure_results:
        print(result)
    
    # Тест импортов
    import_results = test_imports()
    for result in import_results:
        print(result)
    
    # Тест критических проблем
    critical_results = test_critical_issues()
    for result in critical_results:
        print(result)
    
    # Итоги
    print("\n" + "=" * 60)
    print("📊 ИТОГИ:")
    
    all_results = structure_results + import_results + critical_results
    
    ok_count = sum(1 for r in all_results if r.startswith("✅"))
    error_count = sum(1 for r in all_results if r.startswith("❌"))
    warning_count = sum(1 for r in all_results if r.startswith("⚠️"))
    
    print(f"✅ Успешно: {ok_count}")
    print(f"❌ Ошибок: {error_count}")
    print(f"⚠️ Предупреждений: {warning_count}")
    
    if error_count > 0:
        print("\n❗ Рекомендуется запустить автоматические исправления:")
        print("   python diagnostics\\auto_fix.py")
    else:
        print("\n✨ Все основные проверки пройдены!")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
