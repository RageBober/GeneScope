"""
BioForge Complete Diagnostic and Fix Runner
Запускает полную диагностику и предлагает автоматические исправления
"""

import sys
import subprocess
from pathlib import Path

def main():
    print("=" * 70)
    print("🧬 BioForge Project Diagnostic & Repair Tool")
    print("=" * 70)
    
    diagnostics_dir = Path(__file__).parent
    project_root = diagnostics_dir.parent
    
    print(f"\n📁 Проект: {project_root}")
    
    # Шаг 1: Запуск диагностики
    print("\n" + "=" * 70)
    print("ЭТАП 1: ДИАГНОСТИКА")
    print("=" * 70)
    
    diagnostic_script = diagnostics_dir / "full_diagnostic.py"
    print("\n🔍 Запуск полной диагностики...")
    
    try:
        result = subprocess.run(
            [sys.executable, str(diagnostic_script)],
            capture_output=True,
            text=True,
            cwd=str(project_root)
        )
        
        print(result.stdout)
        
        if result.stderr:
            print("Предупреждения:", result.stderr)
            
    except Exception as e:
        print(f"❌ Ошибка при запуске диагностики: {e}")
        return
        
    # Шаг 2: Предложение исправлений
    print("\n" + "=" * 70)
    print("ЭТАП 2: АВТОМАТИЧЕСКИЕ ИСПРАВЛЕНИЯ")
    print("=" * 70)
    
    response = input("\n❓ Применить автоматические исправления? (y/n): ")
    
    if response.lower() == 'y':
        fix_script = diagnostics_dir / "auto_fix.py"
        print("\n🔧 Применение исправлений...")
        
        try:
            result = subprocess.run(
                [sys.executable, str(fix_script)],
                capture_output=True,
                text=True,
                cwd=str(project_root)
            )
            
            print(result.stdout)
            
            if result.stderr:
                print("Предупреждения:", result.stderr)
                
        except Exception as e:
            print(f"❌ Ошибка при применении исправлений: {e}")
            return
            
        print("\n✅ Исправления применены успешно!")
        
    # Шаг 3: Финальные рекомендации
    print("\n" + "=" * 70)
    print("📋 СЛЕДУЮЩИЕ ШАГИ")
    print("=" * 70)
    
    print("""
1. 📦 Установка зависимостей:
   pip install -r requirements.txt
   
2. 🧪 Запуск тестов:
   python -m pytest tests/ -v
   
3. 🚀 Запуск приложения:
   - GUI: python -m genoscope.main --gui
   - API: uvicorn genoscope.api.main:app --reload
   
4. 📚 Документация:
   - Отчет диагностики: diagnostics/report.json
   - Отчет исправлений: diagnostics/fixes_report.txt
   - Бэкап: diagnostics/backup/
   
5. ⚠️ Проверьте:
   - Все изменения перед коммитом
   - Работоспособность критических функций
   - Совместимость с вашей версией Python
""")
    
    print("=" * 70)
    print("✨ Диагностика завершена!")
    print("=" * 70)


if __name__ == "__main__":
    main()
