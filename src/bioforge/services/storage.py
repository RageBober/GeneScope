"""Storage service for file operations."""

import hashlib
import shutil
from pathlib import Path
from typing import BinaryIO
from uuid import uuid4

import aiofiles

from bioforge.config import settings
from bioforge.core.exceptions import StorageError, FileTooLargeError


class StorageService:
    """Service for managing file storage."""
    
    def __init__(self):
        self.base_path = settings.storage.data_dir
        self.samples_path = settings.storage.samples_dir
        self.results_path = settings.storage.results_dir
        self.temp_path = settings.storage.temp_dir
        
        # Ensure directories exist
        for path in [self.base_path, self.samples_path, self.results_path, self.temp_path]:
            path.mkdir(parents=True, exist_ok=True)
    
    def get_sample_path(self, project_id: int, sample_id: int) -> Path:
        """Get storage path for a sample."""
        path = self.samples_path / str(project_id) / str(sample_id)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    def get_results_path(self, project_id: int, sample_id: int, job_id: int) -> Path:
        """Get storage path for job results."""
        path = self.results_path / str(project_id) / str(sample_id) / str(job_id)
        path.mkdir(parents=True, exist_ok=True)
        return path
    
    async def save_file(
        self,
        file: BinaryIO,
        destination: Path,
        max_size: int | None = None,
        chunk_size: int = 1024 * 1024,  # 1MB
    ) -> tuple[Path, int, str]:
        """
        Save uploaded file with size validation.
        
        Returns:
            Tuple of (path, size, md5_hash)
        """
        max_size = max_size or settings.storage.max_file_size
        destination.parent.mkdir(parents=True, exist_ok=True)
        
        hasher = hashlib.md5()
        total_size = 0
        
        async with aiofiles.open(destination, 'wb') as out:
            while True:
                chunk = await file.read(chunk_size)
                if not chunk:
                    break
                
                total_size += len(chunk)
                if total_size > max_size:
                    await out.close()
                    destination.unlink(missing_ok=True)
                    raise FileTooLargeError(total_size, max_size)
                
                hasher.update(chunk)
                await out.write(chunk)
        
        return destination, total_size, hasher.hexdigest()
    
    def generate_unique_filename(self, original_name: str, prefix: str = "") -> str:
        """Generate unique filename preserving extension."""
        suffix = Path(original_name).suffix
        unique_id = uuid4().hex[:8]
        if prefix:
            return f"{prefix}_{unique_id}{suffix}"
        return f"{unique_id}{suffix}"
    
    def delete_file(self, path: Path) -> bool:
        """Delete a file."""
        try:
            if path.exists():
                path.unlink()
                return True
            return False
        except Exception as e:
            raise StorageError(f"Failed to delete file: {e}", str(path))
    
    def delete_directory(self, path: Path) -> bool:
        """Delete a directory and all contents."""
        try:
            if path.exists():
                shutil.rmtree(path)
                return True
            return False
        except Exception as e:
            raise StorageError(f"Failed to delete directory: {e}", str(path))
    
    def get_file_size(self, path: Path) -> int:
        """Get file size in bytes."""
        if not path.exists():
            raise StorageError("File not found", str(path))
        return path.stat().st_size
    
    def file_exists(self, path: Path) -> bool:
        """Check if file exists."""
        return path.exists() and path.is_file()
    
    def list_files(self, directory: Path, pattern: str = "*") -> list[Path]:
        """List files in directory matching pattern."""
        if not directory.exists():
            return []
        return list(directory.glob(pattern))


# Global instance
storage_service = StorageService()
