#!/usr/bin/env python3
"""
Сравнение производительности Poetry vs UV vs pip
"""

import time
import subprocess
import sys
import tempfile
from pathlib import Path
import shutil

def measure_time(func):
    """Декоратор для измерения времени выполнения"""
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        return end - start, result
    return wrapper

@measure_time
def test_pip_install(packages, venv_path):
    """Тест установки через pip"""
    try:
        # Создаем виртуальное окружение
        subprocess.run([sys.executable, "-m", "venv", venv_path], 
                      capture_output=True, check=True)
        
        # Путь к pip в виртуальном окружении
        if sys.platform == "win32":
            pip_path = venv_path / "Scripts" / "pip"
        else:
            pip_path = venv_path / "bin" / "pip"
        
        # Устанавливаем пакеты
        subprocess.run([str(pip_path), "install"] + packages, 
                      capture_output=True, check=True)
        return True
    except Exception as e:
        print(f"pip error: {e}")
        return False

@measure_time
def test_uv_install(packages, venv_path):
    """Тест установки через UV"""
    try:
        # Создаем виртуальное окружение через UV
        subprocess.run(["uv", "venv", str(venv_path)], 
                      capture_output=True, check=True)
        
        # Устанавливаем пакеты
        subprocess.run(["uv", "pip", "install"] + packages, 
                      capture_output=True, check=True,
                      env={**subprocess.os.environ, "VIRTUAL_ENV": str(venv_path)})
        return True
    except Exception as e:
        print(f"UV error: {e}")
        return False

@measure_time
def test_poetry_install(packages, project_path):
    """Тест установки через Poetry"""
    try:
        # Создаем новый Poetry проект
        subprocess.run(["poetry", "new", "test_project", "--no-interaction"], 
                      capture_output=True, check=True, cwd=project_path)
        
        project_dir = project_path / "test_project"
        
        # Добавляем пакеты
        for package in packages:
            subprocess.run(["poetry", "add", package, "--no-interaction"], 
                          capture_output=True, check=True, cwd=project_dir)
        return True
    except Exception as e:
        print(f"Poetry error: {e}")
        return False

def main():
    print("=" * 60)
    print("🏁 СРАВНЕНИЕ ПРОИЗВОДИТЕЛЬНОСТИ: Poetry vs UV vs pip")
    print("=" * 60)
    print()
    
    # Пакеты для тестирования
    test_packages = ["pandas", "numpy", "scikit-learn", "matplotlib"]
    print(f"📦 Тестовые пакеты: {', '.join(test_packages)}")
    print()
    
    results = {}
    
    # Проверяем наличие инструментов
    tools_available = {
        "pip": True,  # pip всегда доступен
        "uv": False,
        "poetry": False
    }
    
    # Проверка UV
    try:
        subprocess.run(["uv", "--version"], capture_output=True, check=True)
        tools_available["uv"] = True
    except:
        print("⚠️ UV не установлен. Установите: pip install uv")
    
    # Проверка Poetry
    try:
        subprocess.run(["poetry", "--version"], capture_output=True, check=True)
        tools_available["poetry"] = True
    except:
        print("⚠️ Poetry не установлен. Установите: pip install poetry")
    
    print()
    print("🔧 Доступные инструменты:")
    for tool, available in tools_available.items():
        status = "✅" if available else "❌"
        print(f"   {status} {tool}")
    
    print()
    print("⏱️ Начинаем тестирование...")
    print("-" * 40)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # Тест pip
        if tools_available["pip"]:
            print("\n📦 Тестирование pip...")
            venv_path = temp_path / "venv_pip"
            time_taken, success = test_pip_install(test_packages, venv_path)
            if success:
                results["pip"] = time_taken
                print(f"   ✅ pip: {time_taken:.2f} секунд")
            else:
                print("   ❌ pip: ошибка установки")
        
        # Тест UV
        if tools_available["uv"]:
            print("\n⚡ Тестирование UV...")
            venv_path = temp_path / "venv_uv"
            time_taken, success = test_uv_install(test_packages, venv_path)
            if success:
                results["uv"] = time_taken
                print(f"   ✅ UV: {time_taken:.2f} секунд")
            else:
                print("   ❌ UV: ошибка установки")
        
        # Тест Poetry
        if tools_available["poetry"]:
            print("\n🎭 Тестирование Poetry...")
            time_taken, success = test_poetry_install(test_packages, temp_path)
            if success:
                results["poetry"] = time_taken
                print(f"   ✅ Poetry: {time_taken:.2f} секунд")
            else:
                print("   ❌ Poetry: ошибка установки")
    
    # Результаты
    print()
    print("=" * 60)
    print("📊 РЕЗУЛЬТАТЫ")
    print("=" * 60)
    
    if not results:
        print("❌ Не удалось провести тесты")
        return
    
    # Сортируем по времени
    sorted_results = sorted(results.items(), key=lambda x: x[1])
    
    print()
    print("🏆 Рейтинг по скорости:")
    for i, (tool, time_taken) in enumerate(sorted_results, 1):
        if i == 1:
            emoji = "🥇"
        elif i == 2:
            emoji = "🥈"
        else:
            emoji = "🥉"
        print(f"{emoji} {i}. {tool.upper()}: {time_taken:.2f} сек")
    
    print()
    print("📈 Сравнение скорости:")
    if "uv" in results and len(results) > 1:
        uv_time = results["uv"]
        for tool, time_taken in results.items():
            if tool != "uv":
                speedup = time_taken / uv_time
                print(f"   UV быстрее {tool} в {speedup:.1f}x раз")
    
    print()
    print("💡 Выводы:")
    if "uv" in results and results["uv"] == min(results.values()):
        print("   ⚡ UV - самый быстрый менеджер пакетов!")
        print("   🚀 Рекомендуется для использования в проекте")
    
    print()
    print("📦 Установка UV:")
    print("   pip install uv")
    print()
    print("🔄 Миграция на UV:")
    print("   python migrate_to_uv.py")
    print()
    print("=" * 60)

if __name__ == "__main__":
    main()
