"""
Модуль валидации данных и файлов для безопасности.
"""
import hashlib
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


class DataValidator:
    """Класс для валидации данных и файлов."""

    ALLOWED_EXTENSIONS = {
        ".csv", ".tsv", ".vcf", ".bam", ".fasta", ".fa", ".fna",
        ".gff", ".gff3", ".gtf", ".json", ".xlsx", ".xls", ".hdf5", ".h5"
    }

    MAX_FILE_SIZE = 500 * 1024 * 1024  # 500MB

    SUSPICIOUS_PATTERNS = [
        "..", "~", "|", ";", "&",
        "<script", "javascript:", "file://", "http://", "https://"
    ]

    @classmethod
    def validate_file_path(cls, file_path: str) -> tuple[bool, str]:
        """Комплексная валидация пути к файлу."""
        try:
            path = Path(file_path).resolve()

            if not path.exists():
                return False, f"Файл не существует: {file_path}"

            if not path.is_file():
                return False, f"Путь не указывает на файл: {file_path}"

            size = path.stat().st_size
            if size > cls.MAX_FILE_SIZE:
                size_mb = size / 1024 / 1024
                return False, f"Файл слишком большой: {size_mb:.1f}MB"

            if size == 0:
                return False, "Файл пустой"

            if path.suffix.lower() not in cls.ALLOWED_EXTENSIONS:
                return False, f"Неподдерживаемое расширение: {path.suffix}"

            path_str = str(path)
            for pattern in cls.SUSPICIOUS_PATTERNS:
                if pattern in path_str:
                    return False, f"Подозрительный паттерн: {pattern}"

            return True, "Валидация прошла успешно"

        except Exception as e:
            return False, f"Ошибка при валидации: {e}"

    @classmethod
    def validate_dataframe(cls, df, min_rows: int = 1, max_rows: int = 1_000_000) -> tuple[bool, str]:
        """Валидация DataFrame."""
        if df is None:
            return False, "DataFrame is None"

        if df.empty:
            return False, "DataFrame пуст"

        if len(df) < min_rows:
            return False, f"Слишком мало строк: {len(df)}"

        if len(df) > max_rows:
            return False, f"Слишком много строк: {len(df)}"

        return True, "DataFrame валиден"

    @classmethod
    def get_file_hash(cls, file_path: str) -> str:
        """Генерация хеша файла для кэширования."""
        path = Path(file_path)
        stat = path.stat()
        content = f"{path.name}:{stat.st_size}:{stat.st_mtime}"
        return hashlib.md5(content.encode()).hexdigest()
