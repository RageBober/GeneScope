"""
Примеры использования GenoScope API
"""
import requests
from pathlib import Path

# Базовый URL API (измените при необходимости)
BASE_URL = "http://localhost:8000"

def test_health():
    """Проверка состояния API"""
    response = requests.get(f"{BASE_URL}/health")
    print("Health check:", response.json())
    return response.status_code == 200

def upload_dataset(file_path: str):
    """Загрузка датасета"""
    with open(file_path, 'rb') as f:
        files = {'file': (Path(file_path).name, f, 'text/csv')}
        response = requests.post(f"{BASE_URL}/datasets/upload", files=files)
    
    if response.status_code == 200:
        data = response.json()
        print(f"Dataset uploaded successfully, ID: {data['dataset_id']}")
        print(f"Rows: {data['rows']}, Columns: {data['cols']}")
        return data['dataset_id']
    else:
        print(f"Upload failed: {response.status_code} - {response.text}")
        return None

def filter_variants(dataset_id: int, **filters):
    """Фильтрация вариантов"""
    params = {"dataset_id": dataset_id}
    params.update(filters)
    
    response = requests.post(f"{BASE_URL}/variants/filter", params=params)
    
    if response.status_code == 200:
        data = response.json()
        print(f"Found {data['total']} variants")
        print(f"Preview: {len(data['preview'])} rows")
        return data
    else:
        print(f"Filter failed: {response.status_code} - {response.text}")
        return None

def main():
    """Основной пример использования"""
    print("🧬 GenoScope API Examples")
    print("=" * 50)
    
    # 1. Проверка здоровья
    if not test_health():
        print("❌ API недоступен")
        return
    
    # 2. Пример создания тестового CSV
    test_csv = Path("test_data.csv")
    if not test_csv.exists():
        csv_content = """CHROM,POS,REF,ALT,GENE,AF,QUAL
1,100,A,T,GENE1,0.01,30
2,200,T,C,GENE2,0.05,40
3,300,G,A,GENE3,0.1,50
X,400,C,G,GENE4,0.5,60
"""
        test_csv.write_text(csv_content)
        print(f"✅ Создан тестовый файл: {test_csv}")
    
    # 3. Загрузка датасета
    dataset_id = upload_dataset(str(test_csv))
    if not dataset_id:
        print("❌ Не удалось загрузить датасет")
        return
    
    # 4. Примеры фильтрации
    print("\n📊 Примеры фильтрации:")
    
    # Фильтр по частоте аллелей
    print("\n1. Редкие варианты (AF < 0.05):")
    filter_variants(dataset_id, max_af=0.05)
    
    # Фильтр по качеству
    print("\n2. Высокое качество (QUAL >= 40):")
    filter_variants(dataset_id, min_qual=40)
    
    # Фильтр по хромосоме
    print("\n3. Только X хромосома:")
    filter_variants(dataset_id, chroms=["X"])
    
    # Комбинированный фильтр с отчетом
    print("\n4. Комбинированный фильтр + отчет:")
    result = filter_variants(
        dataset_id, 
        min_af=0.01, 
        max_af=0.1,
        min_qual=30,
        create_report=True,
        locale="ru",
        clinic="Тестовая клиника"
    )
    
    if result and result.get('report_id'):
        print(f"📄 Отчет создан с ID: {result['report_id']}")
        print(f"🔗 Ссылка: {BASE_URL}/report/{result['report_id']}")
    
    print("\n✅ Все примеры выполнены!")

if __name__ == "__main__":
    main()
