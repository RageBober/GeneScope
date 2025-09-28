# 🧬 GenoScope - Quick Start Guide

## 🚀 Быстрый запуск

GenoScope - это платформа для полного цикла анализа геномных данных, от сырых ридов до аннотированных вариантов.

### 📋 Требования

- Python 3.9+
- 8GB RAM минимум (16GB рекомендуется)
- 50GB свободного места на диске
- Linux/WSL/macOS (Windows через WSL)

### 🎯 Быстрая установка и запуск

#### Вариант 1: Автоматическая установка (рекомендуется)

**Linux/WSL/macOS:**
```bash
# Сделать скрипт исполняемым
chmod +x quick_start.sh

# Запустить полную установку
./quick_start.sh --setup

# Запустить сервер
./quick_start.sh --start
```

**Windows:**
```cmd
# Запустить полную установку
quick_start.bat

# Выбрать опцию 1 для установки
# Затем опцию 2 для запуска сервера
```

**Python (универсальный):**
```bash
# Полная установка
python quick_start.py --setup

# Запуск сервера
python quick_start.py --start
```

#### Вариант 2: Docker (самый простой)

```bash
# Запуск с Docker Compose
docker-compose up -d

# Остановка
docker-compose down
```

#### Вариант 3: Ручная установка

```bash
# 1. Создать виртуальное окружение
python3 -m venv .venv

# 2. Активировать окружение
source .venv/bin/activate  # Linux/macOS
# или
.venv\Scripts\activate.bat  # Windows

# 3. Установить зависимости
pip install -r requirements.txt

# 4. Создать .env файл
cp .env.example .env

# 5. Запустить сервер
uvicorn src.genoscope.api.main:app --reload --host 0.0.0.0 --port 8000
```

### 🌐 Доступ к системе

После запуска GenoScope доступен по адресам:

- **API**: http://localhost:8000
- **Документация**: http://localhost:8000/docs
- **UI**: http://localhost:8000/ui
- **Swagger**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc

### 🧪 Тестовый запуск анализа

```python
# demo_pipeline.py
from pathlib import Path
from src.genoscope.pipeline import PipelineOrchestrator, PipelineConfig

# Настройка pipeline
config = PipelineConfig(
    analysis_type="wgs",
    reference_genome="hg38",
    alignment_tool="bwa",
    variant_caller="gatk",
    threads=8
)

# Запуск анализа
orchestrator = PipelineOrchestrator(config)
result = orchestrator.run_pipeline(
    fastq_r1=Path("data/test/sample_R1.fastq.gz"),
    fastq_r2=Path("data/test/sample_R2.fastq.gz"),
    sample_name="demo_sample"
)

print(f"Статус: {result.status}")
print(f"Найдено вариантов: {result.variant_stats.total_variants}")
```

### 📦 Установка биоинформатических инструментов

GenoScope требует следующие инструменты для полноценной работы:

**Через Conda (рекомендуется):**
```bash
# Установить Miniconda если нет
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Установить инструменты
conda install -c bioconda bwa samtools bcftools fastqc fastp gatk4 minimap2

# Или через скрипт
./quick_start.sh --install-tools
```

**Через apt (Ubuntu/Debian):**
```bash
sudo apt-get update
sudo apt-get install -y bwa samtools bcftools fastqc
```

### 📊 Основные возможности

#### 1. Quality Control
- FastQC анализ
- Адаптер тримминг (fastp)
- Удаление дубликатов
- HTML отчеты

#### 2. Alignment
- BWA для Illumina
- Minimap2 для long-reads
- STAR для RNA-seq
- Автоматическое индексирование

#### 3. Variant Calling
- GATK HaplotypeCaller
- bcftools
- FreeBayes
- Фильтрация вариантов

#### 4. Annotation
- VEP/SnpEff
- ClinVar интеграция
- gnomAD частоты
- Приоритизация вариантов

### 🔧 Конфигурация

Основные настройки в файле `.env`:

```env
# Окружение
GENOSCOPE_ENV=development

# API
API_PORT=8000

# База данных
DATABASE_URL=sqlite:///./genoscope.db

# Redis (для кэша и очередей)
REDIS_URL=redis://localhost:6379/0

# Анализ
DEFAULT_THREADS=8
DEFAULT_MEMORY_GB=16
REFERENCE_GENOME=hg38

# Безопасность (измените в production!)
SECRET_KEY=your-secret-key-change-this
JWT_SECRET_KEY=your-jwt-secret
```

### 📝 API Endpoints

#### Загрузка данных
```bash
# Загрузить FASTQ файл
curl -X POST "http://localhost:8000/datasets/upload" \
  -H "accept: application/json" \
  -H "Content-Type: multipart/form-data" \
  -F "file=@sample.fastq.gz"
```

#### Запуск анализа
```bash
# Запустить полный pipeline
curl -X POST "http://localhost:8000/pipeline/run" \
  -H "Content-Type: application/json" \
  -d '{
    "sample_name": "test_sample",
    "fastq_r1": "sample_R1.fastq.gz",
    "fastq_r2": "sample_R2.fastq.gz",
    "analysis_type": "wgs"
  }'
```

#### Получение результатов
```bash
# Получить статус pipeline
curl "http://localhost:8000/pipeline/status/{pipeline_id}"

# Скачать VCF файл
curl "http://localhost:8000/results/{pipeline_id}/vcf"
```

### 🐛 Решение проблем

#### Ошибка: "Python not found"
```bash
# Установить Python 3.9+
sudo apt-get install python3.9 python3.9-venv python3-pip
```

#### Ошибка: "Port 8000 already in use"
```bash
# Найти процесс
sudo lsof -i :8000

# Остановить процесс
kill -9 <PID>

# Или использовать другой порт
uvicorn src.genoscope.api.main:app --port 8001
```

#### Ошибка: "Module not found"
```bash
# Переустановить зависимости
pip install -r requirements.txt

# Или установить конкретный модуль
pip install fastapi uvicorn pandas biopython
```

#### Ошибка: "BWA/GATK not found"
```bash
# Установить через conda
conda install -c bioconda bwa gatk4

# Или скачать вручную
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xjf bwa-0.7.17.tar.bz2
cd bwa-0.7.17 && make
sudo cp bwa /usr/local/bin/
```

### 📚 Документация

- **API Docs**: http://localhost:8000/docs
- **ReDoc**: http://localhost:8000/redoc
- **GitHub**: [GenoScope Repository](https://github.com/yourusername/genoscope)

### 🤝 Поддержка

- **Email**: support@genoscope.io
- **Telegram**: @genoscope_support
- **Issues**: [GitHub Issues](https://github.com/yourusername/genoscope/issues)

### 📄 Лицензия

MIT License - см. файл LICENSE

---

**GenoScope v1.0** - Современная платформа для геномного анализа 🧬
