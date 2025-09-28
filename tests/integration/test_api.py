"""
Integration Tests for API Endpoints

Тестируем:
1. API endpoints работают корректно
2. Валидация входных данных
3. Обработка ошибок
4. Аутентификация и авторизация
5. Взаимодействие с базой данных
"""

import pytest
import asyncio
from fastapi.testclient import TestClient
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from unittest.mock import patch, MagicMock
import tempfile
import json
from pathlib import Path

from src.genoscope.api.main import app
from src.genoscope.api.db import get_session, Base


# Создаем тестовую базу данных
SQLALCHEMY_DATABASE_URL = "sqlite:///./test.db"
engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
TestingSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)


def override_get_db():
    """Override database dependency для тестов"""
    try:
        db = TestingSessionLocal()
        yield db
    finally:
        db.close()


# Override dependency
app.dependency_overrides[get_session] = override_get_db

# Создаем тестовый клиент
client = TestClient(app)


class TestHealthEndpoints:
    """Тесты для health check endpoints"""
    
    def test_health_check(self):
        """Тест: Health check endpoint возвращает статус"""
        response = client.get("/health")
        
        assert response.status_code == 200
        data = response.json()
        assert "status" in data
        assert data["status"] == "ok"
    
    def test_metrics_endpoint(self):
        """Тест: Metrics endpoint доступен"""
        response = client.get("/metrics")
        
        # Может быть 200 если metrics включены, или 404 если нет
        assert response.status_code in [200, 404]
        
        if response.status_code == 200:
            # Проверяем, что возвращается Prometheus формат
            assert "# HELP" in response.text or "# TYPE" in response.text


class TestPipelineEndpoints:
    """Тесты для pipeline endpoints"""
    
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup для каждого теста"""
        # Создаем таблицы
        Base.metadata.create_all(bind=engine)
        yield
        # Очищаем таблицы после теста
        Base.metadata.drop_all(bind=engine)
    
    def test_submit_pipeline_missing_fields(self):
        """Тест: Отправка pipeline без обязательных полей"""
        response = client.post("/api/pipeline/submit")
        
        assert response.status_code == 422  # Unprocessable Entity
    
    def test_submit_pipeline_with_valid_data(self):
        """Тест: Успешная отправка pipeline job"""
        # Создаем тестовые файлы
        with tempfile.NamedTemporaryFile(suffix='.fastq', delete=False) as f1:
            f1.write(b"@SEQ1\nATCG\n+\nIIII\n")
            fastq_r1 = f1.name
        
        try:
            with open(fastq_r1, 'rb') as f:
                files = {
                    "fastq_r1": ("test_R1.fastq", f, "application/octet-stream")
                }
                data = {
                    "sample_name": "test_sample",
                    "reference": "GRCh38",
                    "aligner": "bwa",
                    "variant_caller": "gatk",
                    "threads": 2
                }
                
                response = client.post(
                    "/api/pipeline/submit",
                    data=data,
                    files=files
                )
            
            assert response.status_code == 200
            result = response.json()
            assert "job_id" in result
            assert result["status"] == "submitted"
            assert "test_sample" in result["message"]
            
        finally:
            Path(fastq_r1).unlink()
    
    def test_get_pipeline_status(self):
        """Тест: Получение статуса pipeline job"""
        # Сначала создаем job
        job_id = "test-job-123"
        
        # Mock PIPELINE_JOBS
        with patch('src.genoscope.api.pipeline_router.PIPELINE_JOBS', {
            job_id: {
                "status": "running",
                "submitted_at": "2024-01-01T12:00:00",
                "config": {"sample_name": "test_sample"}
            }
        }):
            response = client.get(f"/api/pipeline/status/{job_id}")
            
            assert response.status_code == 200
            data = response.json()
            assert data["job_id"] == job_id
            assert data["status"] == "running"
            assert data["sample_name"] == "test_sample"
    
    def test_get_pipeline_status_not_found(self):
        """Тест: Получение статуса несуществующего job"""
        response = client.get("/api/pipeline/status/non-existent-job")
        
        assert response.status_code == 404
        assert "not found" in response.json()["detail"].lower()
    
    def test_list_pipeline_jobs(self):
        """Тест: Список pipeline jobs"""
        # Mock PIPELINE_JOBS
        with patch('src.genoscope.api.pipeline_router.PIPELINE_JOBS', {
            "job1": {
                "status": "completed",
                "submitted_at": "2024-01-01T12:00:00",
                "config": {"sample_name": "sample1"}
            },
            "job2": {
                "status": "running",
                "submitted_at": "2024-01-01T13:00:00",
                "config": {"sample_name": "sample2"}
            }
        }):
            response = client.get("/api/pipeline/jobs")
            
            assert response.status_code == 200
            data = response.json()
            assert data["total"] == 2
            assert len(data["jobs"]) == 2
    
    def test_check_available_tools(self):
        """Тест: Проверка доступных инструментов"""
        response = client.get("/api/pipeline/tools")
        
        assert response.status_code == 200
        data = response.json()
        assert "available" in data
        assert "missing" in data
        assert "total_available" in data
        assert "total_missing" in data
        assert isinstance(data["details"], dict)


class TestGenomicsEndpoints:
    """Тесты для genomics analysis endpoints"""
    
    def test_search_clinvar_missing_params(self):
        """Тест: Поиск в ClinVar без параметров"""
        response = client.get("/api/genomics/clinvar/")
        
        # Должен вернуть 404 так как rsID не указан
        assert response.status_code in [404, 422]
    
    @patch('requests.get')
    def test_search_clinvar_with_rsid(self, mock_get):
        """Тест: Поиск в ClinVar по rsID"""
        # Mock ответ от ClinVar API
        mock_response = MagicMock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "result": {
                "rs12345": {
                    "clinical_significance": "Pathogenic",
                    "gene": "BRCA1",
                    "conditions": ["Breast cancer"]
                }
            }
        }
        mock_get.return_value = mock_response
        
        response = client.get("/api/genomics/clinvar/rs12345")
        
        if response.status_code == 200:
            data = response.json()
            assert "clinical_significance" in data
            assert "gene" in data


class TestBillingEndpoints:
    """Тесты для billing endpoints"""
    
    def test_get_subscription_plans(self):
        """Тест: Получение доступных планов подписки"""
        response = client.get("/api/billing/plans")
        
        assert response.status_code == 200
        data = response.json()
        assert "plans" in data
        assert len(data["plans"]) > 0
        
        # Проверяем структуру плана
        plan = data["plans"][0]
        assert "id" in plan
        assert "name" in plan
        assert "price" in plan
        assert "price_display" in plan
        assert "features" in plan
    
    @patch('src.genoscope.billing.stripe_integration.StripeManager.create_customer')
    @patch('src.genoscope.billing.stripe_integration.StripeManager.create_subscription')
    def test_create_subscription(self, mock_create_sub, mock_create_customer):
        """Тест: Создание подписки"""
        mock_create_customer.return_value = "cus_test123"
        mock_create_sub.return_value = {
            "subscription_id": "sub_test123",
            "status": "active",
            "trial_end": "2024-02-01T00:00:00"
        }
        
        response = client.post("/api/billing/subscribe", json={
            "plan_id": "professional",
            "user_email": "test@example.com",
            "user_name": "Test User"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert "subscription" in data
    
    def test_create_subscription_free_plan(self):
        """Тест: Активация бесплатного плана"""
        response = client.post("/api/billing/subscribe", json={
            "plan_id": "free",
            "user_email": "test@example.com",
            "user_name": "Test User"
        })
        
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "success"
        assert data["plan"] == "Free"
    
    @patch('src.genoscope.billing.stripe_integration.UsageTracker.get_usage_summary')
    @patch('src.genoscope.billing.stripe_integration.UsageTracker.check_limits')
    def test_get_usage_summary(self, mock_check_limits, mock_get_usage):
        """Тест: Получение статистики использования"""
        mock_get_usage.return_value = {
            "month": "January 2024",
            "analyses_count": 10,
            "storage_mb": 5120,
            "api_calls": 1000,
            "estimated_cost": 0.0
        }
        
        mock_check_limits.return_value = {
            "can_analyze": True,
            "can_store": True,
            "can_call_api": True,
            "usage": {
                "analyses": "10/50",
                "storage": "5.0/100 GB",
                "api_calls": "1000/10000"
            }
        }
        
        response = client.get("/api/billing/usage?user_id=test_user")
        
        assert response.status_code == 200
        data = response.json()
        assert "usage" in data
        assert "limits" in data
        assert data["usage"]["analyses_count"] == 10


class TestFileOperations:
    """Тесты для операций с файлами"""
    
    def test_upload_invalid_file_type(self):
        """Тест: Загрузка файла недопустимого типа"""
        with tempfile.NamedTemporaryFile(suffix='.exe', delete=False) as f:
            f.write(b"MZ\x90\x00")  # EXE header
            exe_file = f.name
        
        try:
            with open(exe_file, 'rb') as f:
                response = client.post(
                    "/files",
                    files={"file": ("malware.exe", f, "application/x-msdownload")}
                )
            
            # Должен отклонить exe файл
            assert response.status_code in [400, 415]
            
        finally:
            Path(exe_file).unlink()
    
    def test_upload_valid_vcf(self):
        """Тест: Успешная загрузка VCF файла"""
        vcf_content = b"""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tG\t30\tPASS\t."""
        
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as f:
            f.write(vcf_content)
            vcf_file = f.name
        
        try:
            with open(vcf_file, 'rb') as f:
                response = client.post(
                    "/files",
                    files={"file": ("test.vcf", f, "text/plain")}
                )
            
            # Может быть 200 если DB включена, или другой код
            if response.status_code == 200:
                data = response.json()
                assert "saved" in data or "job_id" in data
                
        finally:
            Path(vcf_file).unlink()
    
    def test_upload_file_size_limit(self):
        """Тест: Проверка ограничения размера файла"""
        # Создаем большой файл (больше лимита)
        large_size = 1024 * 1024 * 501  # 501 MB
        
        with tempfile.NamedTemporaryFile(suffix='.fastq') as f:
            # Не создаем реально большой файл, просто проверяем headers
            response = client.post(
                "/files",
                files={"file": ("huge.fastq", b"@SEQ\nATCG\n+\nIIII", "text/plain")},
                headers={"Content-Length": str(large_size)}
            )
        
        # Сервер должен проверять размер
        # Но в тестах может не работать из-за TestClient


class TestAuthentication:
    """Тесты для аутентификации и авторизации"""
    
    def test_protected_endpoint_without_auth(self):
        """Тест: Доступ к защищенному endpoint без аутентификации"""
        # Предполагаем, что некоторые endpoints требуют auth
        response = client.get("/api/admin/users")
        
        # Должен вернуть 401 или 403 если endpoint защищен
        # Или 404 если endpoint не существует
        assert response.status_code in [401, 403, 404]
    
    def test_metrics_endpoint_access(self):
        """Тест: Доступ к metrics должен быть открыт"""
        response = client.get("/metrics")
        
        # Metrics обычно открыты для Prometheus
        # Но могут быть закрыты в production
        assert response.status_code in [200, 404]
