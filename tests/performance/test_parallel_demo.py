#!/usr/bin/env python3
"""Демонстрация параллельной обработки BioForge с производительными тестами."""

import pandas as pd
import numpy as np
from pathlib import Path
import time
import sys
import argparse

# Добавляем src в путь для импорта
sys.path.insert(0, str(Path(__file__).parent / "src"))

def create_test_genomic_data(size_mb: int = 150, file_type: str = "csv") -> Path:
    """Создает тестовые геномные данные заданного размера.
    
    Args:
        size_mb: Размер файла в МБ
        file_type: Тип файла (csv или vcf)
        
    Returns:
        Путь к созданному файлу
    """
    print(f"📊 Создание тестовых данных {size_mb}MB ({file_type.upper()})...")
    
    # Оценка количества строк (примерно 100 байт на строку)
    rows_needed = int((size_mb * 1024 * 1024) / 100)
    
    if file_type.lower() == "vcf":
        # VCF данные (геномные варианты)
        data = {
            'CHROM': np.random.choice(['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', 'X', 'Y'], rows_needed),
            'POS': np.random.randint(1000000, 250000000, rows_needed),
            'ID': [f'rs{i}' for i in range(rows_needed)],
            'REF': np.random.choice(['A', 'T', 'G', 'C'], rows_needed),
            'ALT': np.random.choice(['A', 'T', 'G', 'C'], rows_needed),
            'QUAL': np.random.uniform(10, 100, rows_needed),
            'FILTER': np.random.choice(['PASS', 'LowQual', 'PASS', 'PASS'], rows_needed),
            'INFO': [f'DP={dp};AF={af:.3f}' for dp, af in zip(
                np.random.randint(10, 200, rows_needed),
                np.random.uniform(0.01, 0.99, rows_needed)
            )],
            'FORMAT': ['GT:DP:GQ'] * rows_needed,
            'SAMPLE1': [f'{gt}:{dp}:{gq}' for gt, dp, gq in zip(
                np.random.choice(['0/0', '0/1', '1/1'], rows_needed),
                np.random.randint(10, 200, rows_needed),
                np.random.randint(20, 99, rows_needed)
            )]
        }
        test_file = Path(f"test_genomic_data_{size_mb}mb.vcf")
        df = pd.DataFrame(data)
        df.to_csv(test_file, sep='\t', index=False)
        
    else:  # CSV
        # CSV данные (общие биологические данные)
        data = {
            'SampleID': [f'Sample_{i:06d}' for i in range(rows_needed)],
            'Gene': np.random.choice([f'Gene_{i}' for i in range(1, 1001)], rows_needed),
            'Expression': np.random.lognormal(3, 1, rows_needed),
            'Chromosome': np.random.choice([f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY'], rows_needed),
            'Position': np.random.randint(1000000, 250000000, rows_needed),
            'QualityScore': np.random.uniform(0, 100, rows_needed),
            'Coverage': np.random.randint(1, 1000, rows_needed),
            'Strand': np.random.choice(['+', '-'], rows_needed),
            'Treatment': np.random.choice(['Control', 'Treated'], rows_needed),
            'BatchEffect': np.random.normal(1, 0.1, rows_needed)
        }
        test_file = Path(f"test_genomic_data_{size_mb}mb.csv")
        df = pd.DataFrame(data)
        df.to_csv(test_file, index=False)
    
    actual_size = test_file.stat().st_size / 1024 / 1024
    print(f"✅ Создан файл: {test_file} ({actual_size:.1f}MB, {len(df)} строк)")
    return test_file


def benchmark_sequential_processing(test_file: Path, file_type: str):
    """Тестирует последовательную обработку."""
    print(f"\n📊 Последовательная обработка ({file_type.upper()}):")
    print("-" * 50)
    
    try:
        from genoscope.main import GenoScopeProcessor
        
        processor = GenoScopeProcessor()
        
        start_time = time.time()
        success = processor.load_file(str(test_file), file_type)
        load_time = time.time() - start_time
        
        if success and processor.data is not None:
            records_count = len(processor.data)
            speed = records_count / load_time if load_time > 0 else 0
            
            print(f"⏱️  Время загрузки: {load_time:.2f} секунд")
            print(f"📊 Загружено записей: {records_count:,}")
            print(f"⚡ Скорость: {speed:.0f} записей/сек")
            
            # Тестируем очистку данных
            start_time = time.time()
            clean_success = processor.clean_data()
            clean_time = time.time() - start_time
            
            if clean_success:
                print(f"🧹 Очистка данных: {clean_time:.2f} секунд")
            
            return {
                'success': True,
                'load_time': load_time,
                'clean_time': clean_time,
                'total_time': load_time + clean_time,
                'records': records_count,
                'speed': speed
            }
        else:
            print("❌ Ошибка при загрузке данных")
            return {'success': False}
            
    except Exception as e:
        print(f"❌ Ошибка последовательной обработки: {e}")
        return {'success': False, 'error': str(e)}


def benchmark_parallel_processing(test_file: Path, file_type: str, workers: int = 8):
    """Тестирует параллельную обработку."""
    print(f"\n⚡ Параллельная обработка ({file_type.upper()}, {workers} воркеров):")
    print("-" * 50)
    
    try:
        from genoscope.parallel import DaskGenomicProcessor
        
        processor = DaskGenomicProcessor(
            n_workers=workers,
            memory_limit="2GB",
            threads_per_worker=2
        )
        
        if not processor.is_distributed:
            print("⚠️  Dask кластер недоступен, используется fallback режим")
        else:
            print(f"🎯 Dask dashboard: {processor.client.dashboard_link}")
        
        start_time = time.time()
        
        # Запускаем параллельную обработку
        results = processor.process_large_file_parallel(
            file_path=test_file,
            file_type=file_type,
            analysis_type="comprehensive",
            chunk_size_mb=50
        )
        
        total_time = time.time() - start_time
        
        if 'error' not in results:
            perf_summary = results.get('performance_summary', {})
            total_records = perf_summary.get('total_records_processed', 0)
            processing_speed = perf_summary.get('average_processing_speed', 0)
            
            print(f"⏱️  Общее время: {total_time:.2f} секунд")
            print(f"📊 Обработано записей: {total_records:,}")
            print(f"⚡ Средняя скорость: {processing_speed:.0f} записей/сек")
            print(f"🔧 Успешные чанки: {results.get('successful_chunks', 0)}/{results.get('total_chunks', 0)}")
            
            # Дополнительные метрики
            if 'performance' in results:
                perf = results['performance']
                print(f"📈 Efficiency score: {perf.get('efficiency_score', 0):.1f}/100")
                print(f"💾 Пиковая память: {perf.get('peak_memory_mb', 0):.1f}MB")
                print(f"🖥️  Среднее CPU: {perf.get('avg_cpu_percent', 0):.1f}%")
            
            processor.close()
            return {
                'success': True,
                'total_time': total_time,
                'records': total_records,
                'speed': processing_speed,
                'chunks': results.get('total_chunks', 0),
                'results': results
            }
        else:
            print(f"❌ Ошибка параллельной обработки: {results['error']}")
            processor.close()
            return {'success': False, 'error': results['error']}
            
    except Exception as e:
        print(f"❌ Ошибка при инициализации параллельной обработки: {e}")
        return {'success': False, 'error': str(e)}


def benchmark_chunk_managers(test_file: Path, file_type: str):
    """Тестирует работу chunk managers."""
    print(f"\n🔧 Тестирование Chunk Managers ({file_type.upper()}):")
    print("-" * 50)
    
    try:
        from genoscope.parallel import get_chunk_manager
        
        chunk_manager = get_chunk_manager(file_type)
        print(f"📦 Используется менеджер: {chunk_manager.__class__.__name__}")
        
        start_time = time.time()
        chunks = chunk_manager.create_chunks(test_file, chunk_size_mb=50)
        chunk_time = time.time() - start_time
        
        total_rows = sum(len(chunk) for chunk in chunks)
        
        print(f"⏱️  Время создания чанков: {chunk_time:.2f} секунд")
        print(f"📊 Создано чанков: {len(chunks)}")
        print(f"📋 Общее количество строк: {total_rows:,}")
        
        # Анализ размеров чанков
        chunk_sizes = [len(chunk) for chunk in chunks]
        print(f"📏 Размеры чанков: min={min(chunk_sizes):,}, max={max(chunk_sizes):,}, avg={sum(chunk_sizes)/len(chunk_sizes):.0f}")
        
        return {
            'success': True,
            'chunk_time': chunk_time,
            'chunks_count': len(chunks),
            'total_rows': total_rows,
            'chunk_sizes': chunk_sizes
        }
        
    except Exception as e:
        print(f"❌ Ошибка тестирования chunk managers: {e}")
        return {'success': False, 'error': str(e)}


def generate_performance_report(seq_results, par_results, chunk_results, test_file: Path):
    """Генерирует итоговый отчет о производительности."""
    print("\n📊 ИТОГОВЫЙ ОТЧЕТ ПРОИЗВОДИТЕЛЬНОСТИ")
    print("=" * 70)
    
    file_size_mb = test_file.stat().st_size / 1024 / 1024
    print(f"📁 Файл: {test_file.name} ({file_size_mb:.1f}MB)")
    
    if seq_results.get('success') and par_results.get('success'):
        seq_time = seq_results['total_time']
        par_time = par_results['total_time']
        speedup = seq_time / par_time if par_time > 0 else 0
        
        print("\n⏱️  ВРЕМЯ ВЫПОЛНЕНИЯ:")
        print(f"   Последовательно: {seq_time:.2f} секунд")
        print(f"   Параллельно:     {par_time:.2f} секунд")
        print(f"   Ускорение:       {speedup:.2f}x")
        
        print("\n⚡ ПРОИЗВОДИТЕЛЬНОСТЬ:")
        print(f"   Последовательно: {seq_results.get('speed', 0):.0f} записей/сек")
        print(f"   Параллельно:     {par_results.get('speed', 0):.0f} записей/сек")
        
        efficiency = "Высокая" if speedup > 3 else "Средняя" if speedup > 1.5 else "Низкая"
        print(f"\n📈 ЭФФЕКТИВНОСТЬ: {efficiency}")
        
        if speedup > 1:
            time_saved = seq_time - par_time
            print(f"💰 ЭКОНОМИЯ ВРЕМЕНИ: {time_saved:.2f} секунд ({(time_saved/seq_time*100):.1f}%)")
    
    if chunk_results.get('success'):
        print("\n🔧 CHUNK MANAGEMENT:")
        print(f"   Время создания чанков: {chunk_results['chunk_time']:.2f} сек")
        print(f"   Количество чанков:     {chunk_results['chunks_count']}")
        print(f"   Средний размер чанка:  {chunk_results['total_rows']/chunk_results['chunks_count']:.0f} строк")


def export_detailed_report(par_results, output_path: str = "performance_report.html"):
    """Экспортирует детальный отчет в HTML."""
    if par_results.get('success') and 'results' in par_results:
        try:
            from genoscope.parallel import PerformanceMonitor
            
            # Создаем mock monitor для демонстрации
            monitor = PerformanceMonitor()
            results = par_results['results']
            
            if 'performance' in results:
                # Создаем HTML отчет
                html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>BioForge Parallel Processing Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        .header {{ background: #667eea; color: white; padding: 20px; border-radius: 8px; }}
        .metric {{ margin: 10px 0; padding: 10px; background: #f8f9fa; border-left: 4px solid #667eea; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 BioForge Parallel Processing Report</h1>
        <p>Generated on {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="metric">
        <h3>📊 Processing Summary</h3>
        <p>Total Time: {par_results['total_time']:.2f} seconds</p>
        <p>Records Processed: {par_results['records']:,}</p>
        <p>Processing Speed: {par_results['speed']:.0f} records/sec</p>
        <p>Chunks Processed: {par_results['chunks']}</p>
    </div>
</body>
</html>
                """
                
                with open(output_path, 'w', encoding='utf-8') as f:
                    f.write(html_content)
                
                print(f"📄 Детальный отчет сохранен: {output_path}")
                
        except Exception as e:
            print(f"⚠️ Не удалось создать детальный отчет: {e}")


def main():
    """Основная функция демонстрации."""
    parser = argparse.ArgumentParser(
        description="Демонстрация параллельной обработки BioForge",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Примеры использования:
  python test_parallel_demo.py                    # Базовый тест CSV (150MB)
  python test_parallel_demo.py --size 300 --vcf   # VCF тест (300MB)
  python test_parallel_demo.py --workers 16       # 16 параллельных воркеров
        """
    )
    
    parser.add_argument("--size", type=int, default=150,
                       help="Размер тестового файла в MB (по умолчанию: 150)")
    parser.add_argument("--vcf", action="store_true",
                       help="Использовать VCF формат вместо CSV")
    parser.add_argument("--workers", type=int, default=8,
                       help="Количество параллельных воркеров (по умолчанию: 8)")
    parser.add_argument("--skip-sequential", action="store_true",
                       help="Пропустить последовательное тестирование")
    parser.add_argument("--output", type=str, default="performance_report.html",
                       help="Файл для детального отчета")
    
    args = parser.parse_args()
    
    print("🧬 BioForge Parallel Processing Demo")
    print("=" * 70)
    
    file_type = "vcf" if args.vcf else "csv"
    
    try:
        # 1. Создание тестовых данных
        test_file = create_test_genomic_data(args.size, file_type)
        
        # 2. Тестирование chunk managers
        chunk_results = benchmark_chunk_managers(test_file, file_type)
        
        # 3. Последовательная обработка (опционально)
        if not args.skip_sequential:
            seq_results = benchmark_sequential_processing(test_file, file_type)
        else:
            seq_results = {'success': False}
        
        # 4. Параллельная обработка
        par_results = benchmark_parallel_processing(test_file, file_type, args.workers)
        
        # 5. Генерация отчета
        generate_performance_report(seq_results, par_results, chunk_results, test_file)
        
        # 6. Экспорт детального отчета
        export_detailed_report(par_results, args.output)
        
        print("\n🎉 Демонстрация завершена!")
        print(f"📁 Тестовый файл: {test_file}")
        print(f"📄 Отчет: {args.output}")
        
        # Cleanup
        cleanup = input("\n🗑️  Удалить тестовый файл? (y/n): ").lower().strip()
        if cleanup == 'y':
            test_file.unlink()
            print(f"🗑️  Файл {test_file} удален")
        
    except KeyboardInterrupt:
        print("\n\n⚠️  Демонстрация прервана пользователем")
    except Exception as e:
        print(f"\n❌ Неожиданная ошибка: {e}")
        import traceback
        traceback.print_exc()
    finally:
        print("\n👋 До свидания!")


if __name__ == "__main__":
    main()
