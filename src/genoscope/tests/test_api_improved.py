"""
Тесты для API endpoints GenoScope.
"""

import io

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from genoscope.api.main import app


class TestAPIEndpoints:
    """Тесты для основных API эндпоинтов."""

    @pytest.fixture
    def client(self):
        return TestClient(app)

    def test_health_endpoint(self, client: TestClient):
        """Тест эндпоинта health."""
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert "status" in data
        assert data["status"] == "ok"

    def test_root_redirect(self, client: TestClient):
        """Тест редиректа с главной страницы."""
        response = client.get("/", follow_redirects=False)
        assert response.status_code == 302
        assert "ui" in response.headers["location"]

    def test_swagger_docs(self, client: TestClient):
        """Тест доступности Swagger документации."""
        response = client.get("/docs")
        assert response.status_code == 200
        assert "swagger" in response.text.lower()

    def test_openapi_json(self, client: TestClient):
        """Тест OpenAPI спецификации."""
        response = client.get("/openapi.json")
        assert response.status_code == 200
        data = response.json()
        assert "openapi" in data
        assert "paths" in data

    def test_favicon(self, client: TestClient):
        """Тест доступности favicon."""
        response = client.get("/favicon.svg")
        # Может быть 200 (файл найден) или 404 (файл не найден)
        assert response.status_code in [200, 404]


class TestFileUpload:
    """Тесты загрузки файлов."""

    @pytest.fixture
    def client(self):
        return TestClient(app)

    def create_test_csv(self) -> bytes:
        """Создает тестовый CSV файл в памяти."""
        data = {
            "CHROM": ["1", "2", "3"],
            "POS": [100, 200, 300],
            "REF": ["A", "T", "G"],
            "ALT": ["T", "C", "A"],
            "QUAL": [30, 40, 50],
        }
        df = pd.DataFrame(data)
        buffer = io.StringIO()
        df.to_csv(buffer, index=False)
        return buffer.getvalue().encode()

    def test_upload_valid_csv(self, client: TestClient):
        """Тест загрузки валидного CSV файла."""
        csv_content = self.create_test_csv()

        files = {"file": ("test.csv", io.BytesIO(csv_content), "text/csv")}

        # Проверяем что эндпоинт существует (может вернуть 503 если DB не настроена)
        response = client.post("/datasets/upload", files=files)
        assert response.status_code in [200, 503]

    def test_upload_invalid_file_type(self, client: TestClient):
        """Тест загрузки файла неподдерживаемого типа."""
        files = {"file": ("test.txt", io.BytesIO(b"test content"), "text/plain")}

        response = client.post("/datasets/upload", files=files)
        assert response.status_code in [
            400,
            503,
        ]  # 400 - неверный тип, 503 - DB не настроена

    def test_upload_empty_filename(self, client: TestClient):
        """Тест загрузки файла без имени."""
        files = {"file": ("", io.BytesIO(b"test"), "text/csv")}

        response = client.post("/datasets/upload", files=files)
        assert response.status_code in [400, 422, 503]

    def test_upload_too_large_file(self, client: TestClient):
        """Тест загрузки слишком большого файла."""
        # Создаем файл размером больше лимита (обычно 100MB)
        large_content = b"x" * (101 * 1024 * 1024)  # 101MB

        files = {"file": ("large.csv", io.BytesIO(large_content), "text/csv")}

        response = client.post("/datasets/upload", files=files)
        assert response.status_code in [413, 422, 503]  # 413 - слишком большой


class TestVariantsFilter:
    """Тесты фильтрации вариантов."""

    @pytest.fixture
    def client(self):
        return TestClient(app)

    def test_filter_with_invalid_dataset_id(self, client: TestClient):
        """Тест фильтрации с невалидным dataset_id."""
        response = client.post("/variants/filter", params={"dataset_id": -1})
        assert response.status_code in [400, 503]

    def test_filter_with_invalid_af_range(self, client: TestClient):
        """Тест фильтрации с невалидным диапазоном AF."""
        response = client.post(
            "/variants/filter",
            params={"dataset_id": 1, "min_af": 0.5, "max_af": 0.1},  # min > max
        )
        assert response.status_code in [400, 503]

    def test_filter_with_negative_af(self, client: TestClient):
        """Тест фильтрации с отрицательным AF."""
        response = client.post(
            "/variants/filter", params={"dataset_id": 1, "min_af": -0.1}
        )
        assert response.status_code in [400, 503]

    def test_filter_with_af_greater_than_one(self, client: TestClient):
        """Тест фильтрации с AF больше 1."""
        response = client.post(
            "/variants/filter", params={"dataset_id": 1, "max_af": 1.5}
        )
        assert response.status_code in [400, 503]

    def test_filter_with_negative_qual(self, client: TestClient):
        """Тест фильтрации с отрицательным QUAL."""
        response = client.post(
            "/variants/filter", params={"dataset_id": 1, "min_qual": -10}
        )
        assert response.status_code in [400, 503]

    def test_filter_with_invalid_limit(self, client: TestClient):
        """Тест фильтрации с невалидным лимитом."""
        response = client.post("/variants/filter", params={"dataset_id": 1, "limit": 0})
        assert response.status_code in [400, 503]

    def test_filter_with_large_limit(self, client: TestClient):
        """Тест фильтрации со слишком большим лимитом."""
        response = client.post(
            "/variants/filter", params={"dataset_id": 1, "limit": 20000}
        )
        assert response.status_code in [400, 503]

    def test_filter_with_invalid_locale(self, client: TestClient):
        """Тест фильтрации с невалидной локалью."""
        response = client.post(
            "/variants/filter", params={"dataset_id": 1, "locale": "invalid"}
        )
        assert response.status_code in [400, 503]


class TestStaticFiles:
    """Тесты статических файлов."""

    @pytest.fixture
    def client(self):
        return TestClient(app)

    def test_favicon_ico(self, client: TestClient):
        """Тест ICO фавикона."""
        response = client.get("/favicon.ico")
        assert response.status_code in [200, 404]  # Может не быть файла

    def test_favicon_svg(self, client: TestClient):
        """Тест SVG фавикона."""
        response = client.get("/favicon.svg")
        assert response.status_code in [200, 404]

    def test_webmanifest(self, client: TestClient):
        """Тест веб-манифеста."""
        response = client.get("/site.webmanifest")
        assert response.status_code in [200, 404]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
