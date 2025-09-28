#!/usr/bin/env python3
"""
BioForge Module Health Check
Проверка работоспособности основных модулей
"""

import sys
from pathlib import Path

# Добавляем src в путь
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root / "src"))

def check_module(module_name, import_path, test_func=None):
    """Проверка отдельного модуля"""
    try:
        exec(f"from {import_path} import *")
        status = "✅ OK"
        
        # Если есть тестовая функция, выполняем её
        if test_func:
            test_func()
            
    except ImportError as e:
        status = f"❌ Import Error: {e}"
    except Exception as e:
        status = f"⚠️ Error: {e}"
    
    print(f"{module_name:30} {status}")
    return "✅" in status


def test_data_ingestion():
    """Тест загрузки данных"""
    from genoscope.data_analysis.data_ingestion import load_data
    # Просто проверяем что функция существует
    assert callable(load_data)


def test_analysis_core():
    """Тест модуля анализа"""
    from genoscope.data_analysis.analysis_core import extract_pca, analyze_data
    assert callable(extract_pca)
    assert callable(analyze_data)


def test_pipeline():
    """Тест pipeline модулей"""
    from genoscope.pipeline.qc import QCMetrics
    from genoscope.pipeline.variant_calling import VariantStats
    
    # Создаем тестовые объекты
    qc = QCMetrics()
    vs = VariantStats()
    
    assert qc.total_reads == 0
    assert vs.total_variants == 0


def test_api():
    """Тест API модуля"""
    from genoscope.api.main import app
    assert app is not None


def main():
    print("=" * 60)
    print("🧬 BioForge Module Health Check")
    print("=" * 60)
    print()
    
    tests = [
        ("Core: Logging", "genoscope.core.logging_config", None),
        ("Core: Validation", "genoscope.core.validation", None),
        ("Core: Security", "genoscope.core.security", None),
        ("Core: Settings", "genoscope.core.settings", None),
        
        ("Data: Ingestion", "genoscope.data_analysis.data_ingestion", test_data_ingestion),
        ("Data: Cleaning", "genoscope.data_analysis.data_cleaning", None),
        ("Data: Filtering", "genoscope.data_analysis.data_filtering", None),
        ("Data: Analysis Core", "genoscope.data_analysis.analysis_core", test_analysis_core),
        ("Data: Visualization", "genoscope.data_analysis.visualization", None),
        
        ("Pipeline: QC", "genoscope.pipeline.qc", None),
        ("Pipeline: Alignment", "genoscope.pipeline.alignment", None),
        ("Pipeline: Variant Calling", "genoscope.pipeline.variant_calling", None),
        ("Pipeline: Main", "genoscope.pipeline.main_pipeline", test_pipeline),
        
        ("API: Main", "genoscope.api.main", test_api),
        ("API: Schemas", "genoscope.api.schemas", None),
        ("API: Services", "genoscope.api.services", None),
        
        ("ML: Imputation", "genoscope.mlmodel.data_cleaning_AI.ml_imputation", None),
        
        ("Main Module", "genoscope.main", None),
    ]
    
    passed = 0
    failed = 0
    
    for module_name, import_path, test_func in tests:
        if check_module(module_name, import_path, test_func):
            passed += 1
        else:
            failed += 1
    
    print()
    print("=" * 60)
    print(f"📊 РЕЗУЛЬТАТЫ:")
    print(f"   ✅ Успешно: {passed}")
    print(f"   ❌ Ошибок: {failed}")
    
    if failed == 0:
        print()
        print("🎉 Все модули работают корректно!")
        print()
        print("🚀 СЛЕДУЮЩИЕ ШАГИ:")
        print("   1. Запустить API: uvicorn genoscope.api.main:app --reload")
        print("   2. Открыть документацию: http://localhost:8000/docs")
        print("   3. Протестировать загрузку файлов через Web UI")
    else:
        print()
        print("⚠️ Обнаружены проблемы. Рекомендуется:")
        print("   1. Запустить: python diagnostics\\final_fix.py")
        print("   2. Установить зависимости: pip install -r requirements.txt")
        print("   3. Проверить снова")
    
    print("=" * 60)
    
    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
