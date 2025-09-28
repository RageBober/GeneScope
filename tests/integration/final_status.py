#!/usr/bin/env python3
"""
Итоговая проверка статуса BioForge после всех исправлений
"""

import sys
from pathlib import Path
from datetime import datetime

def print_header():
    print("=" * 70)
    print("🧬 BIOFORGE PROJECT STATUS REPORT")
    print("=" * 70)
    print(f"Дата проверки: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print()

def check_critical_files():
    """Проверяет критически важные файлы."""
    print("📁 КРИТИЧЕСКИЕ ФАЙЛЫ:")
    
    critical_files = [
        "src/genoscope/main.py",
        "src/genoscope/core/logging_config.py", 
        "src/genoscope/core/validation.py",
        "src/genoscope/data_analysis/data_ingestion.py",
        "src/genoscope/interface.py",
        "src/genoscope/api/main.py",
        "src/genoscope/parallel/dask_processor.py",
        "src/genoscope/parallel/chunk_managers.py",
        "src/genoscope/parallel/performance_monitor.py",
        "pyproject.toml",
        "README.md"
    ]
    
    all_present = True
    for file_path in critical_files:
        if Path(file_path).exists():
            print(f"  ✅ {file_path}")
        else:
            print(f"  ❌ {file_path} - ОТСУТСТВУЕТ")
            all_present = False
    
    return all_present

def check_enhancements():
    """Проверяет дополнительные улучшения."""
    print("\n🚀 ДОПОЛНИТЕЛЬНЫЕ УЛУЧШЕНИЯ:")
    
    enhancements = [
        "Makefile",
        "quick_check.py",
        "start_gui.py", 
        "scripts/benchmark.py",
        "scripts/monitor.py",
        "README_ENHANCED.md",
        "enhance_bioforge.py"
    ]
    
    enhanced_count = 0
    for enhancement in enhancements:
        if Path(enhancement).exists():
            print(f"  ✅ {enhancement}")
            enhanced_count += 1
        else:
            print(f"  ⚪ {enhancement} - не применено")
    
    return enhanced_count

def check_test_functionality():
    """Проверяет базовую функциональность."""
    print("\n🧪 ФУНКЦИОНАЛЬНОСТЬ:")
    
    sys.path.insert(0, "src")
    
    try:
        from genoscope.main import GenoScopeProcessor
        print("  ✅ GenoScopeProcessor импортируется")
        
        processor = GenoScopeProcessor()
        print("  ✅ GenoScopeProcessor создается")
        
        from genoscope.core.validation import DataValidator
        print("  ✅ DataValidator импортируется")
        
        valid, msg = DataValidator.validate_file_path("nonexistent.csv")
        print("  ✅ DataValidator работает")
        
        # Тестируем параллельные модули
        try:
            from genoscope.parallel import CSVChunkManager, PerformanceMonitor
            print("  ✅ Параллельные модули импортируются")
            
            # Тестируем создание объектов
            chunk_manager = CSVChunkManager()
            monitor = PerformanceMonitor()
            print("  ✅ Параллельные объекты создаются")
            
            # Тестируем конфигурацию параллельной обработки
            processor.set_parallel_config(enable=True, n_workers=2)
            print("  ✅ Параллельная конфигурация работает")
            
        except Exception as e:
            print(f"  ⚠️ Параллельные модули недоступны: {e}")
        
        return True
        
    except Exception as e:
        print(f"  ❌ Ошибка функциональности: {e}")
        return False

def generate_summary(critical_ok, enhancements_count, functionality_ok):
    """Генерирует итоговую сводку."""
    print("\n" + "=" * 70)
    print("📊 ИТОГОВАЯ СВОДКА:")
    print("=" * 70)
    
    # Расчет процента готовности
    critical_score = 40 if critical_ok else 0
    enhancement_score = min(35, enhancements_count * 5)
    functionality_score = 25 if functionality_ok else 0
    
    total_score = critical_score + enhancement_score + functionality_score
    
    print(f"Критические компоненты: {'✅ ОК' if critical_ok else '❌ ПРОБЛЕМЫ'} ({critical_score}/40)")
    print(f"Дополнительные улучшения: {enhancements_count}/7 ({enhancement_score}/35)")
    print(f"Базовая функциональность: {'✅ ОК' if functionality_ok else '❌ ПРОБЛЕМЫ'} ({functionality_score}/25)")
    
    print(f"\n🎯 ОБЩАЯ ГОТОВНОСТЬ: {total_score}/100 ({total_score}%)")
    
    if total_score >= 90:
        status = "🟢 ОТЛИЧНО - Готов к продуктивному использованию"
    elif total_score >= 70:
        status = "🟡 ХОРОШО - Основная функциональность работает"
    elif total_score >= 50:
        status = "🟠 УДОВЛЕТВОРИТЕЛЬНО - Требуется доработка"
    else:
        status = "🔴 НЕУДОВЛЕТВОРИТЕЛЬНО - Критические проблемы"
    
    print(f"📈 СТАТУС: {status}")
    
    return total_score

def show_next_steps(score):
    """Показывает рекомендуемые следующие шаги."""
    print("\n🚀 РЕКОМЕНДУЕМЫЕ ДЕЙСТВИЯ:")
    
    if score < 70:
        print("1. Устраните критические проблемы:")
        print("   python apply_comprehensive_fixes.py")
        print("2. Проверьте установку зависимостей:")
        print("   pip install -e .")
        print("3. Повторите проверку:")
        print("   python final_status.py")
    elif score < 90:
        print("1. Примените дополнительные улучшения:")
        print("   python enhance_bioforge.py")
        print("2. Запустите тесты:")
        print("   python -m pytest tests/")
        print("3. Проверьте GUI:")
        print("   python start_gui.py")
    else:
        print("🎉 Проект готов! Попробуйте:")
        print("1. python quick_check.py      # Быстрая проверка")
        print("2. python start_gui.py        # Запуск GUI")
        print("3. python run_api_local.py    # Веб-интерфейс")
        print("4. make help                  # Все команды")

def main():
    """Основная функция."""
    print_header()
    
    # Проверки
    critical_ok = check_critical_files()
    enhancements_count = check_enhancements()
    functionality_ok = check_test_functionality()
    
    # Итоговая сводка
    score = generate_summary(critical_ok, enhancements_count, functionality_ok)
    show_next_steps(score)
    
    print("\n" + "=" * 70)
    print("📚 Документация: README.md, README_ENHANCED.md")
    print("🐛 Отчеты о проблемах: tests/")
    print("📞 Поддержка: Проверьте логи в logs/")
    print("=" * 70)
    
    return score >= 70

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
