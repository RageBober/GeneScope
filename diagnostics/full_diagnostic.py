#!/usr/bin/env python3
"""
BioForge Project Full Diagnostic Script
Анализирует все проблемы проекта и создает отчет
"""

import os
import sys
import json
import ast
import re
from pathlib import Path
from typing import List, Dict, Any
import importlib.util

class ProjectDiagnostic:
    def __init__(self, project_root: Path):
        self.project_root = project_root
        self.src_path = project_root / "src" / "genoscope"
        self.issues = {
            "critical": [],
            "high": [],
            "medium": [],
            "low": []
        }
        
    def run_full_diagnostic(self):
        """Запуск полной диагностики"""
        print("🔍 Запуск полной диагностики проекта BioForge...")
        print("=" * 60)
        
        # 1. Проверка структуры проекта
        self.check_project_structure()
        
        # 2. Проверка импортов
        self.check_imports()
        
        # 3. Проверка синтаксиса типов
        self.check_type_annotations()
        
        # 4. Проверка логирования
        self.check_logging_config()
        
        # 5. Проверка обработки ошибок
        self.check_error_handling()
        
        # 6. Проверка зависимостей
        self.check_dependencies()
        
        # 7. Проверка GUI компонентов
        self.check_gui_components()
        
        # 8. Проверка анализа данных
        self.check_data_analysis_modules()
        
        # 9. Проверка валидации файлов
        self.check_file_validation()
        
        # 10. Проверка тестов
        self.check_tests()
        
        # Генерация отчета
        self.generate_report()
        
    def check_project_structure(self):
        """Проверка структуры проекта"""
        print("\n📁 Проверка структуры проекта...")
        
        required_dirs = [
            "src/genoscope",
            "src/genoscope/core",
            "src/genoscope/data_analysis",
            "src/genoscope/api",
            "src/genoscope/pipeline",
            "tests",
            "frontend"
        ]
        
        missing_dirs = []
        for dir_path in required_dirs:
            full_path = self.project_root / dir_path
            if not full_path.exists():
                missing_dirs.append(dir_path)
                
        if missing_dirs:
            self.issues["critical"].append({
                "type": "MISSING_DIRECTORIES",
                "details": f"Отсутствуют директории: {missing_dirs}",
                "fix": "Создать недостающие директории"
            })
            
        # Проверка __init__.py файлов
        for root, dirs, files in os.walk(self.src_path):
            if any(f.endswith('.py') and f != '__init__.py' for f in files):
                init_file = Path(root) / '__init__.py'
                if not init_file.exists():
                    rel_path = Path(root).relative_to(self.project_root)
                    self.issues["medium"].append({
                        "type": "MISSING_INIT",
                        "details": f"Отсутствует __init__.py в {rel_path}",
                        "fix": f"Создать пустой __init__.py файл"
                    })
                    
    def check_imports(self):
        """Проверка импортов и модулей"""
        print("\n🔗 Проверка импортов...")
        
        # Проверка отсутствующего interface.py
        interface_path = self.src_path / "interface.py"
        if not interface_path.exists():
            self.issues["critical"].append({
                "type": "MISSING_MODULE",
                "details": "Отсутствует модуль interface.py (GUI на Tkinter)",
                "fix": "Создать interface.py с классом GenoScopeApp или удалить импорты"
            })
            
        # Проверка импортов в main.py
        main_path = self.src_path / "main.py"
        if main_path.exists():
            content = main_path.read_text()
            
            # Проверка импортов из несуществующих модулей
            import_pattern = r'from genoscope\.(\w+) import (\w+)'
            imports = re.findall(import_pattern, content)
            
            for module, item in imports:
                module_path = self.src_path / f"{module}.py"
                if not module_path.exists() and not (self.src_path / module).is_dir():
                    self.issues["high"].append({
                        "type": "BROKEN_IMPORT",
                        "details": f"Импорт из несуществующего модуля: genoscope.{module}",
                        "fix": f"Создать модуль {module} или исправить импорт"
                    })
                    
    def check_type_annotations(self):
        """Проверка совместимости типов"""
        print("\n🏷️ Проверка аннотаций типов...")
        
        py_files = list(self.src_path.rglob("*.py"))
        modern_syntax_count = 0
        
        for py_file in py_files:
            try:
                content = py_file.read_text(encoding='utf-8')
                
                # Проверка современного синтаксиса типов (Python 3.10+)
                if re.search(r'\w+\s*:\s*\w+\s*\|\s*None', content):
                    modern_syntax_count += 1
                    rel_path = py_file.relative_to(self.project_root)
                    self.issues["high"].append({
                        "type": "MODERN_TYPE_SYNTAX",
                        "file": str(rel_path),
                        "details": f"Использование | для Union (Python 3.10+)",
                        "fix": "Заменить на typing.Union или Optional"
                    })
                    
                # Проверка dict[...], list[...] вместо Dict[...], List[...]
                if re.search(r':\s*(dict|list|tuple|set)\[', content):
                    rel_path = py_file.relative_to(self.project_root)
                    self.issues["medium"].append({
                        "type": "LOWERCASE_GENERICS",
                        "file": str(rel_path),
                        "details": "Использование dict[]/list[] вместо Dict[]/List[]",
                        "fix": "Импортировать из typing и использовать капитализированные версии"
                    })
                    
            except Exception as e:
                continue
                
    def check_logging_config(self):
        """Проверка конфигурации логирования"""
        print("\n📝 Проверка логирования...")
        
        logging_config_path = self.src_path / "core" / "logging_config.py"
        if logging_config_path.exists():
            content = logging_config_path.read_text()
            
            # Проверка дублирования handlers
            if '"handlers": ["console", "console_error"]' in content:
                self.issues["high"].append({
                    "type": "DUPLICATE_LOGGING",
                    "details": "Дублирование вывода логов (console и console_error)",
                    "fix": "Использовать фильтры для разделения уровней логирования"
                })
                
            # Проверка отсутствия фильтров
            if 'InfoAndBelowFilter' not in content:
                self.issues["medium"].append({
                    "type": "MISSING_LOG_FILTER",
                    "details": "Отсутствует фильтр для разделения stdout/stderr",
                    "fix": "Добавить InfoAndBelowFilter класс"
                })
                
    def check_error_handling(self):
        """Проверка обработки ошибок"""
        print("\n⚠️ Проверка обработки ошибок...")
        
        py_files = list(self.src_path.rglob("*.py"))
        
        for py_file in py_files:
            try:
                content = py_file.read_text(encoding='utf-8')
                rel_path = py_file.relative_to(self.project_root)
                
                # Проверка на raise SystemExit в библиотечных функциях
                if 'raise SystemExit' in content and 'main.py' not in str(py_file):
                    self.issues["high"].append({
                        "type": "SYSTEMEXIT_IN_LIBRARY",
                        "file": str(rel_path),
                        "details": "SystemExit в библиотечной функции",
                        "fix": "Заменить на ValueError или специфичное исключение"
                    })
                    
                # Проверка на голые except:
                if re.search(r'except\s*:', content):
                    self.issues["medium"].append({
                        "type": "BARE_EXCEPT",
                        "file": str(rel_path),
                        "details": "Использование except без указания типа исключения",
                        "fix": "Указать конкретный тип исключения"
                    })
                    
                # Проверка на pass в except блоках
                if re.search(r'except.*:\s*pass', content):
                    self.issues["low"].append({
                        "type": "SILENT_EXCEPTION",
                        "file": str(rel_path),
                        "details": "Игнорирование исключений без логирования",
                        "fix": "Добавить логирование или обработку"
                    })
                    
            except Exception:
                continue
                
    def check_dependencies(self):
        """Проверка зависимостей"""
        print("\n📦 Проверка зависимостей...")
        
        # Проверка соответствия requirements.txt и pyproject.toml
        req_path = self.project_root / "requirements.txt"
        pyproject_path = self.project_root / "pyproject.toml"
        
        if req_path.exists():
            req_content = req_path.read_text()
            
            # Проверка критических зависимостей для GUI
            if 'tkinter' not in req_content.lower():
                self.issues["medium"].append({
                    "type": "MISSING_DEPENDENCY",
                    "details": "tkinter не указан в зависимостях (нужен для GUI)",
                    "fix": "tkinter встроен в Python, но нужно документировать"
                })
                
    def check_gui_components(self):
        """Проверка GUI компонентов"""
        print("\n🖼️ Проверка GUI компонентов...")
        
        # Проверка наличия interface.py или альтернативных GUI модулей
        gui_modules = [
            self.src_path / "interface.py",
            self.src_path / "gui" / "main.py",
            self.src_path / "ui" / "app.py"
        ]
        
        if not any(module.exists() for module in gui_modules):
            self.issues["critical"].append({
                "type": "MISSING_GUI",
                "details": "Отсутствуют модули GUI (interface.py)",
                "fix": "Создать GUI модуль или использовать только API/CLI"
            })
            
    def check_data_analysis_modules(self):
        """Проверка модулей анализа данных"""
        print("\n📊 Проверка модулей анализа данных...")
        
        analysis_core = self.src_path / "data_analysis" / "analysis_core.py"
        if analysis_core.exists():
            content = analysis_core.read_text()
            
            # Проверка функции PCA
            if 'def extract_pca' in content:
                # Проверка на неправильную валидацию размерности
                if 'X.shape[0] < n_components or X.shape[1] < n_components' in content:
                    self.issues["high"].append({
                        "type": "PCA_VALIDATION_ERROR",
                        "details": "Неправильная проверка размерности в PCA",
                        "fix": "Проверять только X.shape[1] < n_components"
                    })
                    
                # Проверка на отсутствие нормализации
                if 'StandardScaler' not in content:
                    self.issues["medium"].append({
                        "type": "PCA_NO_SCALING",
                        "details": "Отсутствует нормализация данных перед PCA",
                        "fix": "Добавить StandardScaler перед PCA"
                    })
                    
        # Проверка visualization.py
        viz_path = self.src_path / "data_analysis" / "visualization.py"
        if viz_path.exists():
            content = viz_path.read_text()
            
            # Проверка на блокирующие plt.show()
            if 'plt.show()' in content and 'show_plot' not in content:
                self.issues["high"].append({
                    "type": "BLOCKING_PLOT_SHOW",
                    "details": "plt.show() блокирует GUI",
                    "fix": "Добавить параметр show_plot для условного отображения"
                })
                
    def check_file_validation(self):
        """Проверка валидации файлов"""
        print("\n🔒 Проверка валидации файлов...")
        
        # Проверка наличия валидации в data_ingestion.py
        ingestion_path = self.src_path / "data_analysis" / "data_ingestion.py"
        if ingestion_path.exists():
            content = ingestion_path.read_text()
            
            # Проверка на отсутствие валидации размера
            if 'MAX_FILE_SIZE' not in content:
                self.issues["critical"].append({
                    "type": "NO_FILE_SIZE_VALIDATION",
                    "details": "Отсутствует проверка размера загружаемых файлов",
                    "fix": "Добавить MAX_FILE_SIZE и проверку размера"
                })
                
            # Проверка на отсутствие валидации расширений
            if 'ALLOWED_EXTENSIONS' not in content:
                self.issues["critical"].append({
                    "type": "NO_EXTENSION_VALIDATION",
                    "details": "Отсутствует проверка расширений файлов",
                    "fix": "Добавить список разрешенных расширений"
                })
                
            # Проверка на отсутствие проверки путей
            if 'resolve()' not in content and 'expanduser()' not in content:
                self.issues["high"].append({
                    "type": "NO_PATH_RESOLUTION",
                    "details": "Отсутствует разрешение путей (уязвимость path traversal)",
                    "fix": "Использовать Path.resolve() для безопасности"
                })
                
    def check_tests(self):
        """Проверка тестов"""
        print("\n🧪 Проверка тестов...")
        
        tests_dir = self.project_root / "tests"
        if tests_dir.exists():
            test_files = list(tests_dir.glob("test_*.py"))
            
            if len(test_files) == 0:
                self.issues["high"].append({
                    "type": "NO_TESTS",
                    "details": "Отсутствуют тестовые файлы",
                    "fix": "Создать тесты для критической функциональности"
                })
            else:
                # Проверка на placeholder тесты
                for test_file in test_files:
                    content = test_file.read_text()
                    if 'assert True' in content and len(content) < 500:
                        rel_path = test_file.relative_to(self.project_root)
                        self.issues["medium"].append({
                            "type": "PLACEHOLDER_TEST",
                            "file": str(rel_path),
                            "details": "Тест-заглушка без реальной проверки",
                            "fix": "Написать реальные тесты"
                        })
                        
    def generate_report(self):
        """Генерация отчета о найденных проблемах"""
        print("\n" + "=" * 60)
        print("📋 ОТЧЕТ О НАЙДЕННЫХ ПРОБЛЕМАХ")
        print("=" * 60)
        
        total_issues = sum(len(issues) for issues in self.issues.values())
        
        print(f"\n📊 Всего найдено проблем: {total_issues}")
        print(f"   🔴 Критических: {len(self.issues['critical'])}")
        print(f"   🟠 Высокий приоритет: {len(self.issues['high'])}")
        print(f"   🟡 Средний приоритет: {len(self.issues['medium'])}")
        print(f"   🟢 Низкий приоритет: {len(self.issues['low'])}")
        
        # Детальный отчет по категориям
        for priority, emoji in [
            ("critical", "🔴"),
            ("high", "🟠"),
            ("medium", "🟡"),
            ("low", "🟢")
        ]:
            if self.issues[priority]:
                print(f"\n{emoji} {priority.upper()} ПРИОРИТЕТ:")
                print("-" * 50)
                for i, issue in enumerate(self.issues[priority], 1):
                    print(f"\n{i}. [{issue['type']}]")
                    if 'file' in issue:
                        print(f"   📁 Файл: {issue['file']}")
                    print(f"   ❌ Проблема: {issue['details']}")
                    print(f"   ✅ Решение: {issue['fix']}")
                    
        # Сохранение отчета в JSON
        report_path = self.project_root / "diagnostics" / "report.json"
        with open(report_path, 'w', encoding='utf-8') as f:
            json.dump(self.issues, f, ensure_ascii=False, indent=2)
        print(f"\n💾 Отчет сохранен в: {report_path}")
        
        # Рекомендации
        print("\n" + "=" * 60)
        print("🎯 РЕКОМЕНДАЦИИ ПО ИСПРАВЛЕНИЮ")
        print("=" * 60)
        
        print("\n📌 ПЕРВООЧЕРЕДНЫЕ ДЕЙСТВИЯ:")
        print("1. Создать отсутствующий interface.py или удалить импорты")
        print("2. Исправить современный синтаксис типов на совместимый")
        print("3. Добавить валидацию файлов (размер, расширения, пути)")
        print("4. Исправить дублирование в логировании")
        print("5. Заменить SystemExit на ValueError в библиотечных функциях")
        
        print("\n📌 ВАЖНЫЕ ИСПРАВЛЕНИЯ:")
        print("1. Исправить проверку размерности в PCA")
        print("2. Добавить параметр show_plot в визуализации")
        print("3. Добавить нормализацию данных перед PCA")
        print("4. Написать реальные тесты вместо заглушек")
        
        print("\n📌 УЛУЧШЕНИЯ:")
        print("1. Добавить __init__.py где отсутствуют")
        print("2. Улучшить обработку исключений")
        print("3. Добавить документацию")
        print("4. Настроить CI/CD pipeline")


if __name__ == "__main__":
    # Определяем корневую директорию проекта
    project_root = Path(__file__).parent.parent
    
    # Запускаем диагностику
    diagnostic = ProjectDiagnostic(project_root)
    diagnostic.run_full_diagnostic()
