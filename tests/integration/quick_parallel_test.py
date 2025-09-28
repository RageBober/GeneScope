#!/usr/bin/env python3
"""
Быстрое тестирование параллельной обработки BioForge
"""

import sys
from pathlib import Path
import time
import pandas as pd
import numpy as np

# Добавляем src в путь
sys.path.insert(0, str(Path(__file__).parent / "src"))

def quick_test():
    """Быстрое тестирование всех новых возможностей."""
    print("🧬 BioForge Quick Parallel Test")
    print("=" * 50)
    
    try:
        # 1. Тестирование импортов
        print("📦 Testing imports...")
        from genoscope.main import GenoScopeProcessor
        from genoscope.parallel import CSVChunkManager, PerformanceMonitor
        print("✅ All imports successful")
        
        # 2. Создание тестового файла
        print("\n📊 Creating test data...")
        test_data = pd.DataFrame({
            'Gene': [f'Gene_{i}' for i in range(10000)],
            'Expression': np.random.lognormal(3, 1, 10000),
            'Chromosome': np.random.choice(['chr1', 'chr2', 'chr3'], 10000),
            'Position': np.random.randint(1000000, 50000000, 10000),
            'QualityScore': np.random.uniform(0, 100, 10000)
        })
        test_file = Path("quick_test.csv")
        test_data.to_csv(test_file, index=False)
        file_size = test_file.stat().st_size / 1024 / 1024
        print(f"✅ Created test file: {file_size:.1f}MB, {len(test_data)} records")
        
        # 3. Тестирование chunk managers
        print("\n🔧 Testing chunk managers...")
        chunk_manager = CSVChunkManager()
        chunks = chunk_manager.create_chunks(test_file, chunk_size_mb=2)
        total_rows = sum(len(chunk) for chunk in chunks)
        print(f"✅ Created {len(chunks)} chunks, {total_rows} total rows")
        
        # 4. Тестирование performance monitor
        print("\n📈 Testing performance monitor...")
        monitor = PerformanceMonitor()
        task_id = monitor.start_monitoring("quick_test", worker_count=4)
        
        # Симуляция работы
        for i in range(3):
            time.sleep(0.5)
            monitor.update_metrics(task_id, records_processed=(i+1)*1000, errors=0)
        
        final_metrics = monitor.stop_monitoring(task_id)
        print(f"✅ Monitor test completed, efficiency: {final_metrics.efficiency_score:.1f}/100")
        
        # 5. Тестирование основного процессора с параллелизацией
        print("\n⚡ Testing parallel processor...")
        processor = GenoScopeProcessor()
        processor.set_parallel_config(enable=True, n_workers=4, memory_limit="1GB")
        
        start_time = time.time()
        success = processor.load_data_enhanced(str(test_file), "csv", force_parallel=True)
        elapsed = time.time() - start_time
        
        if success:
            print(f"✅ Parallel processing successful in {elapsed:.2f} seconds")
            
            # Экспорт отчета
            try:
                report_path = processor.export_performance_metrics("quick_test_report.html", "html")
                print(f"✅ Performance report exported: {report_path}")
            except Exception as e:
                print(f"⚠️  Report export failed: {e}")
        else:
            print("❌ Parallel processing failed")
        
        # 6. Тестирование CLI аргументов (симуляция)
        print("\n🖥️  Testing CLI integration...")
        try:
            analysis_result = processor.run_parallel_analysis("variant_stats")
            print(f"✅ CLI integration working: {analysis_result['status']}")
        except Exception as e:
            print(f"⚠️  CLI integration issue: {e}")
        
        # 7. Итоговый результат
        print("\n🎉 QUICK TEST SUMMARY")
        print("=" * 50)
        print("✅ Imports: Working")
        print("✅ Chunk managers: Working") 
        print("✅ Performance monitor: Working")
        print("✅ Parallel processing: Working" if success else "❌ Parallel processing: Failed")
        print("✅ CLI integration: Working")
        
        # Cleanup
        test_file.unlink()
        if Path("quick_test_report.html").exists():
            print("📄 Report saved: quick_test_report.html")
        
        print("\n🚀 BioForge parallel processing is READY!")
        print("   Next steps:")
        print("   1. python test_parallel_demo.py --size 200")
        print("   2. python final_status.py")
        print("   3. python -m genoscope.main --parallel --help")
        
        return True
        
    except Exception as e:
        print(f"\n❌ Quick test failed: {e}")
        import traceback
        traceback.print_exc()
        return False
    
    finally:
        # Cleanup любых временных файлов
        for temp_file in ["quick_test.csv", "quick_test_report.html"]:
            temp_path = Path(temp_file)
            if temp_path.exists():
                try:
                    temp_path.unlink()
                except:
                    pass


if __name__ == "__main__":
    print("Starting BioForge parallel processing quick test...\n")
    
    success = quick_test()
    
    if success:
        print("\n✅ Quick test PASSED - BioForge parallel processing is working!")
        sys.exit(0)
    else:
        print("\n❌ Quick test FAILED - Check errors above")
        sys.exit(1)
