"""
Performance Tests using Locust

Тестируем:
1. Нагрузочное тестирование API
2. Проверка масштабируемости
3. Поиск узких мест производительности
4. Стресс-тестирование
"""

from locust import HttpUser, task, between, events
from locust.env import Environment
from locust.stats import stats_printer, stats_history
from locust.log import setup_logging
import random
import json
import time
import gevent


class GenoScopeUser(HttpUser):
    """
    Симулируем поведение пользователя GenoScope
    """
    
    # Пользователи делают запросы каждые 1-5 секунд
    wait_time = between(1, 5)
    
    def on_start(self):
        """Выполняется при старте каждого пользователя"""
        # Симулируем логин
        self.client.post("/api/auth/login", json={
            "email": f"user{random.randint(1, 1000)}@example.com",
            "password": "password123"
        })
        
        # Сохраняем данные пользователя
        self.user_id = f"user_{random.randint(1, 10000)}"
        self.sample_names = [f"sample_{i}" for i in range(5)]
    
    @task(10)
    def view_dashboard(self):
        """Частая операция - просмотр dashboard"""
        with self.client.get("/api/dashboard", 
                            catch_response=True,
                            name="Dashboard") as response:
            if response.status_code == 200:
                response.success()
            else:
                response.failure(f"Got status code {response.status_code}")
    
    @task(5)
    def search_clinvar(self):
        """Поиск в ClinVar - средняя частота"""
        rsids = ["rs80357906", "rs121913343", "rs1799966", "rs12345", "rs67890"]
        rsid = random.choice(rsids)
        
        with self.client.get(f"/api/genomics/clinvar/{rsid}",
                            catch_response=True,
                            name="ClinVar Search") as response:
            if response.status_code in [200, 404]:  # 404 тоже ок для несуществующих
                response.success()
            else:
                response.failure(f"Unexpected status {response.status_code}")
    
    @task(3)
    def upload_vcf(self):
        """Загрузка VCF файла - менее частая операция"""
        # Генерируем небольшой VCF для теста
        vcf_content = self._generate_vcf(num_variants=100)
        
        files = {"file": ("test.vcf", vcf_content, "text/plain")}
        
        with self.client.post("/api/genomics/upload/vcf",
                             files=files,
                             catch_response=True,
                             name="VCF Upload") as response:
            if response.status_code == 200:
                response.success()
                # Сохраняем ID для последующих операций
                data = response.json()
                if "file_id" in data:
                    self.uploaded_files = getattr(self, 'uploaded_files', [])
                    self.uploaded_files.append(data["file_id"])
            else:
                response.failure(f"Upload failed: {response.status_code}")
    
    @task(2)
    def submit_pipeline(self):
        """Отправка pipeline job - редкая но тяжелая операция"""
        data = {
            "sample_name": random.choice(self.sample_names),
            "reference": "GRCh38",
            "aligner": random.choice(["bwa", "minimap2"]),
            "variant_caller": random.choice(["gatk", "freebayes"]),
            "threads": 4
        }
        
        # Создаем mock FASTQ данные
        fastq_content = self._generate_fastq(num_reads=1000)
        files = {
            "fastq_r1": ("sample_R1.fastq", fastq_content, "text/plain"),
            "fastq_r2": ("sample_R2.fastq", fastq_content, "text/plain")
        }
        
        with self.client.post("/api/pipeline/submit",
                             data=data,
                             files=files,
                             catch_response=True,
                             name="Pipeline Submit",
                             timeout=30) as response:
            if response.status_code == 200:
                response.success()
                data = response.json()
                if "job_id" in data:
                    self.job_ids = getattr(self, 'job_ids', [])
                    self.job_ids.append(data["job_id"])
            else:
                response.failure(f"Pipeline submission failed: {response.status_code}")
    
    @task(8)
    def check_job_status(self):
        """Проверка статуса job - частая операция"""
        if hasattr(self, 'job_ids') and self.job_ids:
            job_id = random.choice(self.job_ids)
            
            with self.client.get(f"/api/pipeline/status/{job_id}",
                                catch_response=True,
                                name="Job Status") as response:
                if response.status_code in [200, 404]:
                    response.success()
                else:
                    response.failure(f"Status check failed: {response.status_code}")
    
    @task(4)
    def get_usage_stats(self):
        """Получение статистики использования"""
        with self.client.get(f"/api/billing/usage?user_id={self.user_id}",
                            catch_response=True,
                            name="Usage Stats") as response:
            if response.status_code == 200:
                response.success()
            else:
                response.failure(f"Failed to get usage: {response.status_code}")
    
    @task(1)
    def download_report(self):
        """Скачивание отчета - тяжелая операция"""
        if hasattr(self, 'job_ids') and self.job_ids:
            job_id = random.choice(self.job_ids)
            
            with self.client.get(f"/api/pipeline/download/{job_id}/full_report",
                                catch_response=True,
                                name="Download Report",
                                timeout=60) as response:
                if response.status_code in [200, 404]:
                    response.success()
                else:
                    response.failure(f"Download failed: {response.status_code}")
    
    def _generate_vcf(self, num_variants=100):
        """Генерация тестового VCF файла"""
        vcf = ["##fileformat=VCFv4.2"]
        vcf.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        
        for i in range(num_variants):
            chrom = random.randint(1, 22)
            pos = random.randint(1000000, 50000000)
            ref = random.choice(['A', 'T', 'C', 'G'])
            alt = random.choice(['A', 'T', 'C', 'G'])
            while alt == ref:
                alt = random.choice(['A', 'T', 'C', 'G'])
            qual = random.randint(20, 100)
            
            vcf.append(f"chr{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\tDP=30")
        
        return "\n".join(vcf)
    
    def _generate_fastq(self, num_reads=1000):
        """Генерация тестового FASTQ файла"""
        fastq = []
        bases = ['A', 'T', 'C', 'G']
        
        for i in range(num_reads):
            seq = ''.join(random.choices(bases, k=150))
            qual = 'I' * 150  # Высокое качество
            fastq.extend([
                f"@READ_{i}",
                seq,
                "+",
                qual
            ])
        
        return "\n".join(fastq)


class AdminUser(HttpUser):
    """
    Симулируем администратора - более тяжелые операции
    """
    
    wait_time = between(2, 10)
    weight = 1  # Меньше админов чем обычных пользователей
    
    def on_start(self):
        """Admin login"""
        self.client.post("/api/auth/login", json={
            "email": "admin@genoscope.com",
            "password": "admin_password"
        })
    
    @task(5)
    def view_all_jobs(self):
        """Просмотр всех jobs - тяжелая операция"""
        with self.client.get("/api/pipeline/jobs?limit=100",
                            catch_response=True,
                            name="Admin: All Jobs") as response:
            if response.status_code == 200:
                response.success()
            else:
                response.failure(f"Failed: {response.status_code}")
    
    @task(3)
    def generate_reports(self):
        """Генерация отчетов"""
        report_types = ["usage", "billing", "performance", "errors"]
        
        data = {
            "report_type": random.choice(report_types),
            "start_date": "2024-01-01",
            "end_date": "2024-01-31"
        }
        
        with self.client.post("/api/admin/reports/generate",
                             json=data,
                             catch_response=True,
                             name="Admin: Generate Report",
                             timeout=30) as response:
            if response.status_code in [200, 201]:
                response.success()
            else:
                response.failure(f"Report generation failed: {response.status_code}")
    
    @task(2)
    def manage_users(self):
        """Управление пользователями"""
        with self.client.get("/api/admin/users?page=1&limit=50",
                            catch_response=True,
                            name="Admin: List Users") as response:
            if response.status_code == 200:
                response.success()
            else:
                response.failure(f"Failed: {response.status_code}")
    
    @task(1)
    def system_metrics(self):
        """Получение системных метрик"""
        with self.client.get("/metrics",
                            catch_response=True,
                            name="Admin: System Metrics") as response:
            if response.status_code == 200:
                response.success()
            else:
                response.failure(f"Metrics failed: {response.status_code}")


class StressTestUser(HttpUser):
    """
    Пользователь для стресс-тестирования - максимальная нагрузка
    """
    
    wait_time = between(0.1, 0.5)  # Очень частые запросы
    weight = 0  # По умолчанию отключен, включается только для стресс-тестов
    
    @task
    def hammer_api(self):
        """Постоянные запросы к API"""
        endpoints = [
            "/api/health",
            "/api/genomics/clinvar/rs12345",
            "/api/pipeline/jobs",
            "/api/billing/plans"
        ]
        
        endpoint = random.choice(endpoints)
        
        with self.client.get(endpoint,
                            catch_response=True,
                            name=f"Stress: {endpoint}") as response:
            if response.status_code < 500:
                response.success()
            else:
                response.failure(f"Server error: {response.status_code}")


# Event handlers для сбора метрик
@events.test_start.add_listener
def on_test_start(environment, **kwargs):
    """Начало теста"""
    print("Starting performance test...")
    print(f"Target host: {environment.host}")
    print(f"Total users: {environment.parsed_options.num_users}")
    print(f"Spawn rate: {environment.parsed_options.spawn_rate}")


@events.request.add_listener
def on_request(request_type, name, response_time, response_length, exception, **kwargs):
    """Логирование каждого запроса для анализа"""
    if response_time > 2000:  # Запросы дольше 2 секунд
        print(f"SLOW REQUEST: {name} took {response_time}ms")
    
    if exception:
        print(f"REQUEST FAILED: {name} - {exception}")


@events.test_stop.add_listener
def on_test_stop(environment, **kwargs):
    """Конец теста - выводим итоговую статистику"""
    print("\n" + "="*50)
    print("Performance Test Results")
    print("="*50)
    
    # Общая статистика
    stats = environment.stats
    print(f"Total requests: {stats.total.num_requests}")
    print(f"Total failures: {stats.total.num_failures}")
    print(f"Average response time: {stats.total.avg_response_time:.2f}ms")
    print(f"Min response time: {stats.total.min_response_time:.2f}ms")
    print(f"Max response time: {stats.total.max_response_time:.2f}ms")
    
    # Проверка SLA
    if stats.total.avg_response_time > 500:
        print("⚠️  WARNING: Average response time exceeds 500ms SLA")
    
    if stats.total.num_failures / stats.total.num_requests > 0.01:
        print("⚠️  WARNING: Error rate exceeds 1% SLA")
    
    print("\nTop slowest endpoints:")
    sorted_stats = sorted(stats.entries.values(), 
                         key=lambda x: x.avg_response_time, 
                         reverse=True)[:5]
    
    for stat in sorted_stats:
        print(f"  {stat.name}: {stat.avg_response_time:.2f}ms")


# Функция для запуска программного теста
def run_performance_test(host="http://localhost:8000", users=100, spawn_rate=10, duration=300):
    """
    Запуск performance теста программно
    
    Args:
        host: URL тестируемого сервера
        users: Количество пользователей
        spawn_rate: Скорость создания пользователей (users/sec)
        duration: Длительность теста в секундах
    """
    
    setup_logging("INFO", None)
    
    # Создаем environment
    env = Environment(user_classes=[GenoScopeUser, AdminUser], host=host)
    
    # Запускаем spawn
    env.runner.start(users, spawn_rate=spawn_rate)
    
    # Запускаем статистику
    gevent.spawn(stats_printer(env.stats))
    
    # Ждем
    gevent.spawn_later(duration, lambda: env.runner.quit())
    
    # Запускаем greenlets
    env.runner.greenlet.join()
    
    # Возвращаем статистику
    return {
        "total_requests": env.stats.total.num_requests,
        "total_failures": env.stats.total.num_failures,
        "avg_response_time": env.stats.total.avg_response_time,
        "rps": env.stats.total.current_rps
    }


if __name__ == "__main__":
    # Запуск через CLI: locust -f locustfile.py --host=http://localhost:8000
    # Или программно:
    results = run_performance_test(
        host="http://localhost:8000",
        users=50,
        spawn_rate=5,
        duration=60
    )
    print(f"Test completed: {results}")
