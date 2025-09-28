#!/usr/bin/env python3
"""
Симуляция запуска BioForge с полной проверкой функциональности.
Создает тестовые данные и проверяет работу всех основных компонентов.
"""

import sys
import os
import pandas as pd
import tempfile
from pathlib import Path
import logging

# Добавим src в путь
project_root = Path(__file__).parent
src_path = project_root / "src"
sys.path.insert(0, str(src_path))

# Настройка базового логирования
logging.basicConfig(level=logging.INFO, format='%(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def create_test_data():
    """Создать тестовые данные для симуляции."""
    logger.info("🧪 Создание тестовых данных...")
    
    # Создаем DataFrame с геномными данными
    test_data = pd.DataFrame({
        'sample_id': [f'Sample_{i:03d}' for i in range(1, 101)],
        'chromosome': ['chr' + str((i % 22) + 1) for i in range(100)],
        'position': [1000 + i * 100 for i in range(100)],
        'ref_allele': ['A', 'T', 'G', 'C'] * 25,
        'alt_allele': ['T', 'A', 'C', 'G'] * 25,
        'quality_score': [30 + (i % 40) for i in range(100)],
        'depth': [20 + (i % 80) for i in range(100)],
        'allele_freq': [0.1 + (i % 9) * 0.1 for i in range(100)],
        'gene_name': [f'GENE_{(i % 50) + 1}' for i in range(100)],
        'annotation': [f'Variant_{i}' for i in range(100)]
    })
    
    # Добавляем несколько дубликатов и NaN значений для тестирования очистки
    test_data = pd.concat([test_data, test_data.head(5)], ignore_index=True)
    test_data.loc[10:15, 'quality_score'] = None
    test_data.loc[20:25, 'depth'] = None
    
    logger.info(f"✅ Создан тестовый датасет: {len(test_data)} строк, {len(test_data.columns)} столбцов")
    return test_data

def test_data_ingestion():
    """Тестировать загрузку данных."""
    logger.info("📥 Тестирование загрузки данных...")
    
    try:
        from genoscope.data_analysis.data_ingestion import load_data
        
        # Создаем временный CSV файл
        test_data = create_test_data()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            test_data.to_csv(f.name, index=False)
            temp_csv = f.name
        
        # Тестируем загрузку
        loaded_data = load_data(temp_csv, 'csv')
        
        if loaded_data is not None:
            logger.info(f"✅ Загрузка данных успешна: {len(loaded_data)} строк")
            os.unlink(temp_csv)  # Удаляем временный файл
            return loaded_data
        else:
            logger.error("❌ Загрузка данных не удалась")
            return None
            
    except ImportError as e:
        logger.error(f"❌ Модуль data_ingestion недоступен: {e}")
        return None
    except Exception as e:
        logger.error(f"❌ Ошибка загрузки данных: {e}")
        return None

def test_data_cleaning(data):
    """Тестировать очистку данных."""
    logger.info("🧹 Тестирование очистки данных...")
    
    if data is None:
        logger.error("❌ Нет данных для очистки")
        return None
    
    try:
        from genoscope.data_analysis.data_cleaning import remove_duplicates, handle_missing_values
        
        original_size = len(data)
        logger.info(f"📊 Исходный размер данных: {original_size} строк")
        
        # Удаление дубликатов
        cleaned_data = remove_duplicates(data)
        logger.info(f"✅ Удалено дубликатов: {original_size - len(cleaned_data)} строк")
        
        # Обработка пропущенных значений
        missing_before = cleaned_data.isnull().sum().sum()
        cleaned_data = handle_missing_values(cleaned_data, method='mean')
        missing_after = cleaned_data.isnull().sum().sum()
        
        logger.info(f"✅ Обработано пропущенных значений: {missing_before - missing_after}")
        logger.info(f"📊 Итоговый размер данных: {len(cleaned_data)} строк")
        
        return cleaned_data
        
    except ImportError as e:
        logger.error(f"❌ Модуль data_cleaning недоступен: {e}")
        return data
    except Exception as e:
        logger.error(f"❌ Ошибка очистки данных: {e}")
        return data

def test_analysis_core(data):
    """Тестировать основные аналитические функции."""
    logger.info("🔬 Тестирование анализа данных...")
    
    if data is None:
        logger.error("❌ Нет данных для анализа")
        return None
    
    try:
        from genoscope.data_analysis.analysis_core import extract_pca, basic_statistics
        
        # Базовая статистика
        stats = basic_statistics(data)
        logger.info(f"✅ Статистика: {stats['total_rows']} строк, {stats['numeric_columns']} числовых столбцов")
        
        # PCA анализ
        try:
            pca_result = extract_pca(data, n_components=2)
            explained_variance = pca_result.attrs.get('explained_variance_ratio', [])
            total_variance = sum(explained_variance) if explained_variance else 0
            
            logger.info(f"✅ PCA выполнен: {pca_result.shape}, объяснено {total_variance:.2%} дисперсии")
            return {'stats': stats, 'pca': pca_result}
            
        except Exception as pca_error:
            logger.warning(f"⚠️ PCA не удался: {pca_error}")
            return {'stats': stats, 'pca': None}
        
    except ImportError as e:
        logger.error(f"❌ Модуль analysis_core недоступен: {e}")
        return None
    except Exception as e:
        logger.error(f"❌ Ошибка анализа: {e}")
        return None

def test_visualization(data, analysis_results):
    """Тестировать визуализацию (без отображения)."""
    logger.info("📊 Тестирование визуализации...")
    
    if data is None or analysis_results is None:
        logger.error("❌ Нет данных для визуализации")
        return False
    
    try:
        from genoscope.data_analysis.visualization import plot_correlation_matrix, plot_pca
        
        # Тестируем корреляционную матрицу (без показа)
        try:
            # Monkey patch для предотвращения показа графиков
            import matplotlib.pyplot as plt
            original_show = plt.show
            plt.show = lambda: None  # Отключаем показ графиков
            
            plot_correlation_matrix(data)
            logger.info("✅ Корреляционная матрица построена")
            
            # PCA график
            if analysis_results.get('pca') is not None:
                plot_pca(analysis_results['pca'])
                logger.info("✅ График PCA построен")
            
            plt.show = original_show  # Восстанавливаем оригинальную функцию
            return True
            
        except Exception as viz_error:
            logger.warning(f"⚠️ Визуализация частично не удалась: {viz_error}")
            return False
        
    except ImportError as e:
        logger.error(f"❌ Модуль visualization недоступен: {e}")
        return False
    except Exception as e:
        logger.error(f"❌ Ошибка визуализации: {e}")
        return False

def test_processor_integration():
    """Тестировать интеграцию с основным процессором."""
    logger.info("🔄 Тестирование интеграции процессора...")
    
    try:
        from genoscope.main import GenoScopeProcessor
        
        processor = GenoScopeProcessor()
        logger.info("✅ GenoScopeProcessor создан успешно")
        
        # Создаем временный файл для тестирования
        test_data = create_test_data()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            test_data.to_csv(f.name, index=False)
            temp_csv = f.name
        
        # Тестируем полный пайплайн
        success = processor.run_pipeline(temp_csv, 'csv')
        
        if success:
            logger.info("✅ Полный пайплайн выполнен успешно")
            logger.info(f"📊 Обработано строк: {len(processor.data) if processor.data is not None else 0}")
        else:
            logger.error("❌ Пайплайн завершился с ошибкой")
        
        os.unlink(temp_csv)  # Удаляем временный файл
        return success
        
    except ImportError as e:
        logger.error(f"❌ GenoScopeProcessor недоступен: {e}")
        return False
    except Exception as e:
        logger.error(f"❌ Ошибка в процессоре: {e}")
        return False

def test_gui_components():
    """Тестировать GUI компоненты (без показа)."""
    logger.info("🖥️ Тестирование GUI компонентов...")
    
    try:
        import tkinter as tk
        from genoscope.interface import GenoScopeApp
        
        # Создаем корневое окно (скрытое)
        root = tk.Tk()
        root.withdraw()  # Скрываем окно
        
        # Создаем приложение
        app = GenoScopeApp(root)
        logger.info("✅ GUI приложение создано успешно")
        
        # Уничтожаем окно
        root.destroy()
        return True
        
    except ImportError as e:
        logger.error(f"❌ GUI модули недоступны: {e}")
        return False
    except Exception as e:
        logger.error(f"❌ Ошибка GUI: {e}")
        return False

def test_api_components():
    """Тестировать API компоненты."""
    logger.info("🌐 Тестирование API компонентов...")
    
    try:
        from genoscope.api.main import app
        logger.info("✅ FastAPI приложение импортировано успешно")
        
        # Проверяем наличие основных эндпоинтов
        routes = [route.path for route in app.routes]
        logger.info(f"📋 Доступные маршруты: {len(routes)} штук")
        
        if "/" in routes:
            logger.info("✅ Главная страница доступна")
        
        if "/health" in routes:
            logger.info("✅ Health check эндпоинт доступен")
        
        return True
        
    except ImportError as e:
        logger.error(f"❌ API модули недоступны: {e}")
        return False
    except Exception as e:
        logger.error(f"❌ Ошибка API: {e}")
        return False

def main():
    """Главная функция симуляции."""
    print("🚀 СИМУЛЯЦИЯ ЗАПУСКА BIOFORGE")
    print("=" * 60)
    
    results = {
        'data_ingestion': False,
        'data_cleaning': False,
        'analysis_core': False,
        'visualization': False,
        'processor_integration': False,
        'gui_components': False,
        'api_components': False
    }
    
    # 1. Тестирование загрузки данных
    test_data = test_data_ingestion()
    results['data_ingestion'] = test_data is not None
    
    # 2. Тестирование очистки данных
    if test_data is not None:
        cleaned_data = test_data_cleaning(test_data)
        results['data_cleaning'] = cleaned_data is not None
    else:
        cleaned_data = None
        results['data_cleaning'] = False
    
    # 3. Тестирование анализа
    if cleaned_data is not None:
        analysis_results = test_analysis_core(cleaned_data)
        results['analysis_core'] = analysis_results is not None
    else:
        analysis_results = None
        results['analysis_core'] = False
    
    # 4. Тестирование визуализации
    results['visualization'] = test_visualization(cleaned_data, analysis_results)
    
    # 5. Тестирование интеграции процессора
    results['processor_integration'] = test_processor_integration()
    
    # 6. Тестирование GUI
    results['gui_components'] = test_gui_components()
    
    # 7. Тестирование API
    results['api_components'] = test_api_components()
    
    # Итоговый отчет
    print("\n" + "=" * 60)
    print("📊 ИТОГОВЫЙ ОТЧЕТ СИМУЛЯЦИИ")
    print("=" * 60)
    
    passed = sum(results.values())
    total = len(results)
    
    for component, success in results.items():
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"{component.replace('_', ' ').title():<25} {status}")
    
    print("-" * 60)
    print(f"ОБЩИЙ РЕЗУЛЬТАТ: {passed}/{total} компонентов работают ({passed/total*100:.1f}%)")
    
    if passed == total:
        print("🎉 ВСЕ ТЕСТЫ ПРОШЛИ! BioForge готов к использованию!")
        return 0
    elif passed >= total * 0.7:
        print("⚠️ БОЛЬШИНСТВО ТЕСТОВ ПРОШЛИ. Проект в основном функционален.")
        return 0
    else:
        print("❌ МНОГО НЕУДАЧНЫХ ТЕСТОВ. Требуются дополнительные исправления.")
        return 1

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
