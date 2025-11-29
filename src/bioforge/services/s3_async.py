"""
Async S3/MinIO Storage Service.

Uses aioboto3 for async operations with S3-compatible storage.
"""

from contextlib import asynccontextmanager
from io import BytesIO
from pathlib import Path
from typing import Any, AsyncIterator, BinaryIO

import aioboto3
from botocore.config import Config

from bioforge.config import settings
from bioforge.core.logging import get_logger

logger = get_logger(__name__)


class AsyncS3Service:
    """Async service for S3/MinIO object storage operations."""
    
    def __init__(self):
        self._session = aioboto3.Session()
        self._config = Config(
            signature_version='s3v4',
            s3={'addressing_style': 'path'},
            retries={'max_attempts': 3, 'mode': 'adaptive'},
        )
        self._endpoint_url = f"{'https' if settings.minio.secure else 'http'}://{settings.minio.endpoint}"
    
    @asynccontextmanager
    async def _get_client(self):
        """Get S3 client context manager."""
        async with self._session.client(
            's3',
            endpoint_url=self._endpoint_url,
            aws_access_key_id=settings.minio.access_key,
            aws_secret_access_key=settings.minio.secret_key,
            config=self._config,
        ) as client:
            yield client
    
    async def ensure_buckets(self) -> None:
        """Create required buckets if they don't exist."""
        buckets = [
            settings.minio.bucket_raw,
            settings.minio.bucket_aligned,
            settings.minio.bucket_variants,
            settings.minio.bucket_reports,
            settings.minio.bucket_references,
        ]
        
        async with self._get_client() as client:
            # Get existing buckets
            response = await client.list_buckets()
            existing = {b['Name'] for b in response.get('Buckets', [])}
            
            for bucket in buckets:
                if bucket not in existing:
                    try:
                        await client.create_bucket(Bucket=bucket)
                        logger.info("Created bucket", bucket=bucket)
                    except Exception as e:
                        logger.error("Failed to create bucket", bucket=bucket, error=str(e))
    
    async def upload_file(
        self,
        bucket: str,
        key: str,
        file_path: str | Path,
        content_type: str | None = None,
        metadata: dict[str, str] | None = None,
    ) -> str:
        """
        Upload a file to S3.
        
        Returns:
            S3 URI in format s3://bucket/key
        """
        extra_args: dict[str, Any] = {}
        if content_type:
            extra_args['ContentType'] = content_type
        if metadata:
            extra_args['Metadata'] = metadata
        
        async with self._get_client() as client:
            await client.upload_file(
                str(file_path),
                bucket,
                key,
                ExtraArgs=extra_args or None,
            )
        
        logger.info("Uploaded file", bucket=bucket, key=key, path=str(file_path))
        return f"s3://{bucket}/{key}"
    
    async def upload_data(
        self,
        bucket: str,
        key: str,
        data: bytes | BinaryIO,
        content_type: str = "application/octet-stream",
        metadata: dict[str, str] | None = None,
    ) -> str:
        """Upload binary data to S3."""
        if isinstance(data, bytes):
            data = BytesIO(data)
        
        async with self._get_client() as client:
            await client.upload_fileobj(
                data,
                bucket,
                key,
                ExtraArgs={
                    'ContentType': content_type,
                    **(metadata or {}),
                },
            )
        
        logger.info("Uploaded data", bucket=bucket, key=key)
        return f"s3://{bucket}/{key}"
    
    async def download_file(
        self,
        bucket: str,
        key: str,
        file_path: str | Path,
    ) -> Path:
        """Download file from S3 to local path."""
        file_path = Path(file_path)
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        async with self._get_client() as client:
            await client.download_file(bucket, key, str(file_path))
        
        logger.info("Downloaded file", bucket=bucket, key=key, path=str(file_path))
        return file_path
    
    async def get_object(self, bucket: str, key: str) -> bytes:
        """Get object data as bytes."""
        async with self._get_client() as client:
            response = await client.get_object(Bucket=bucket, Key=key)
            async with response['Body'] as stream:
                return await stream.read()
    
    async def delete_object(self, bucket: str, key: str) -> bool:
        """Delete object from S3."""
        try:
            async with self._get_client() as client:
                await client.delete_object(Bucket=bucket, Key=key)
            logger.info("Deleted object", bucket=bucket, key=key)
            return True
        except Exception as e:
            logger.error("Delete failed", bucket=bucket, key=key, error=str(e))
            return False
    
    async def object_exists(self, bucket: str, key: str) -> bool:
        """Check if object exists."""
        try:
            async with self._get_client() as client:
                await client.head_object(Bucket=bucket, Key=key)
            return True
        except Exception:
            return False
    
    async def get_object_info(self, bucket: str, key: str) -> dict[str, Any]:
        """Get object metadata."""
        async with self._get_client() as client:
            response = await client.head_object(Bucket=bucket, Key=key)
            return {
                "bucket": bucket,
                "key": key,
                "size": response.get('ContentLength'),
                "etag": response.get('ETag', '').strip('"'),
                "content_type": response.get('ContentType'),
                "last_modified": response.get('LastModified'),
                "metadata": response.get('Metadata', {}),
            }
    
    async def list_objects(
        self,
        bucket: str,
        prefix: str = "",
        max_keys: int = 1000,
    ) -> list[dict[str, Any]]:
        """List objects in bucket with optional prefix."""
        objects = []
        
        async with self._get_client() as client:
            paginator = client.get_paginator('list_objects_v2')
            
            async for page in paginator.paginate(
                Bucket=bucket,
                Prefix=prefix,
                PaginationConfig={'MaxItems': max_keys},
            ):
                for obj in page.get('Contents', []):
                    objects.append({
                        "key": obj['Key'],
                        "size": obj['Size'],
                        "last_modified": obj['LastModified'],
                        "etag": obj['ETag'].strip('"'),
                    })
        
        return objects
    
    async def generate_presigned_url(
        self,
        bucket: str,
        key: str,
        expires_in: int = 3600,
        method: str = "get_object",
    ) -> str:
        """Generate presigned URL for object access."""
        async with self._get_client() as client:
            return await client.generate_presigned_url(
                method,
                Params={'Bucket': bucket, 'Key': key},
                ExpiresIn=expires_in,
            )
    
    # ─────────────────────────────────────────────────────────────────────
    # Convenience methods for genomic data
    # ─────────────────────────────────────────────────────────────────────
    
    def _sample_key(self, project_id: int, sample_id: int, filename: str) -> str:
        """Generate key for sample file."""
        return f"{project_id}/{sample_id}/{filename}"
    
    async def upload_fastq(
        self,
        project_id: int,
        sample_id: int,
        file_path: Path,
    ) -> str:
        """Upload FASTQ file to raw bucket."""
        key = self._sample_key(project_id, sample_id, file_path.name)
        return await self.upload_file(
            settings.minio.bucket_raw,
            key,
            file_path,
            content_type="application/gzip" if file_path.suffix == ".gz" else "text/plain",
        )
    
    async def upload_bam(
        self,
        project_id: int,
        sample_id: int,
        file_path: Path,
    ) -> str:
        """Upload BAM/CRAM file to aligned bucket."""
        key = self._sample_key(project_id, sample_id, file_path.name)
        return await self.upload_file(
            settings.minio.bucket_aligned,
            key,
            file_path,
            content_type="application/octet-stream",
        )
    
    async def upload_vcf(
        self,
        project_id: int,
        sample_id: int,
        file_path: Path,
    ) -> str:
        """Upload VCF file to variants bucket."""
        key = self._sample_key(project_id, sample_id, file_path.name)
        return await self.upload_file(
            settings.minio.bucket_variants,
            key,
            file_path,
            content_type="application/gzip" if file_path.suffix == ".gz" else "text/plain",
        )
    
    async def get_sample_files(
        self,
        project_id: int,
        sample_id: int,
    ) -> dict[str, list[dict]]:
        """Get all files for a sample across all buckets."""
        prefix = f"{project_id}/{sample_id}/"
        
        return {
            "raw": await self.list_objects(settings.minio.bucket_raw, prefix),
            "aligned": await self.list_objects(settings.minio.bucket_aligned, prefix),
            "variants": await self.list_objects(settings.minio.bucket_variants, prefix),
            "reports": await self.list_objects(settings.minio.bucket_reports, prefix),
        }


# Global instance (lazy init)
_async_s3: AsyncS3Service | None = None


def get_async_s3_service() -> AsyncS3Service:
    """Get async S3 service instance."""
    global _async_s3
    if _async_s3 is None:
        _async_s3 = AsyncS3Service()
    return _async_s3


async def init_s3() -> None:
    """Initialize S3 service and create buckets."""
    service = get_async_s3_service()
    await service.ensure_buckets()
