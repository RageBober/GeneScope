# genoscope/core/cache.py
import pandas as pd
import hashlib
import pickle
from pathlib import Path


class DataCache:
    def __init__(self, cache_dir: Path = Path(".genoscope_cache")):
        self.cache_dir = cache_dir
        self.cache_dir.mkdir(exist_ok=True)

    def get_file_hash(self, file_path: str) -> str:
        """Generate hash for file content and metadata."""
        path = Path(file_path)
        stat = path.stat()
        content = f"{path.name}:{stat.st_size}:{stat.st_mtime}"
        return hashlib.md5(content.encode()).hexdigest()

    def get_cached_data(self, file_hash: str) -> pd.DataFrame | None:
        """Retrieve cached data if exists."""
        cache_file = self.cache_dir / f"{file_hash}.pkl"
        if cache_file.exists():
            try:
                with open(cache_file, "rb") as f:
                    return pickle.load(f)
            except Exception:
                cache_file.unlink()  # Remove corrupted cache
        return None

    def cache_data(self, file_hash: str, data: pd.DataFrame):
        """Cache processed data."""
        cache_file = self.cache_dir / f"{file_hash}.pkl"
        with open(cache_file, "wb") as f:
            pickle.dump(data, f)
