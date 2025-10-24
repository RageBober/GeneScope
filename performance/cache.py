"""
💾 Intelligent Caching System

Redis-based caching for intermediate results:
- Cache expensive computation results
- Automatic expiration
- Memory-efficient serialization
"""

from __future__ import annotations

import hashlib
import logging
import pickle
from typing import Any, Callable, Optional

logger = logging.getLogger(__name__)

# Try importing Redis
try:
    import redis
    REDIS_AVAILABLE = True
except ImportError:
    REDIS_AVAILABLE = False
    logger.warning("Redis not available. Install with: pip install redis")


class CacheManager:
    """
    Intelligent cache manager using Redis.

    Examples:
        >>> cache = CacheManager()
        >>> result = cache.get_or_compute(
        ...     "expensive_operation",
        ...     lambda: expensive_function(),
        ...     ttl=3600
        ... )
    """

    def __init__(
        self,
        host: str = "localhost",
        port: int = 6379,
        db: int = 0,
        ttl: int = 3600,
    ):
        """
        Initialize cache manager.

        Args:
            host: Redis host
            port: Redis port
            db: Redis database number
            ttl: Default time-to-live (seconds)
        """
        self.ttl = ttl
        self.redis_available = REDIS_AVAILABLE

        if self.redis_available:
            try:
                self.client = redis.Redis(
                    host=host,
                    port=port,
                    db=db,
                    decode_responses=False
                )
                self.client.ping()
                logger.info(f"✅ Redis connected: {host}:{port}")
            except Exception as exc:
                logger.error(f"Redis connection failed: {exc}")
                self.redis_available = False

    def get_or_compute(
        self,
        key: str,
        compute_func: Callable[[], Any],
        ttl: Optional[int] = None,
    ) -> Any:
        """
        Get from cache or compute if not exists.

        Args:
            key: Cache key
            compute_func: Function to compute value if not cached
            ttl: Time-to-live (seconds)

        Returns:
            Cached or computed value
        """
        if not self.redis_available:
            return compute_func()

        ttl = ttl or self.ttl

        # Try to get from cache
        cached = self.get(key)
        if cached is not None:
            logger.info(f"✅ Cache hit: {key}")
            return cached

        # Compute
        logger.info(f"⚠️  Cache miss: {key}, computing...")
        value = compute_func()

        # Store in cache
        self.set(key, value, ttl=ttl)

        return value

    def get(self, key: str) -> Optional[Any]:
        """Get value from cache"""
        if not self.redis_available:
            return None

        try:
            data = self.client.get(key)
            if data:
                return pickle.loads(data)
        except Exception as exc:
            logger.error(f"Cache get error: {exc}")

        return None

    def set(self, key: str, value: Any, ttl: Optional[int] = None):
        """Set value in cache"""
        if not self.redis_available:
            return

        ttl = ttl or self.ttl

        try:
            data = pickle.dumps(value)
            self.client.setex(key, ttl, data)
            logger.info(f"💾 Cached: {key} (TTL={ttl}s)")
        except Exception as exc:
            logger.error(f"Cache set error: {exc}")

    def clear(self):
        """Clear all cache"""
        if self.redis_available:
            self.client.flushdb()
            logger.info("Cache cleared")
