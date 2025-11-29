"""S3/MinIO storage service."""

import logging
from io import BytesIO
from pathlib import Path
from typing import BinaryIO

from minio import Minio
from minio.error import S3Error

from bioforge.config import settings

logger = logging.getLogger("bioforge.storage.s3")


class S3StorageService:
    """Service for S3/MinIO object storage operations."""
    
    def __init__(self):
        self.client = Minio(
            settings.minio.endpoint,
            access_key=settings.minio.access_key,
            secret_key=settings.minio.secret_key,
            secure=settings.minio.secure,
        )
        self._ensure_buckets()
    
    def _ensure_buckets(self) -> None:
        """Create required buckets if they don't exist."""
        buckets = [
            settings.minio.bucket_raw,
            settings.minio.bucket_aligned,
            settings.minio.bucket_variants,
            settings.minio.bucket_reports,
            settings.minio.bucket_references,
        ]
        
        for bucket in buckets:
            try:
                if not self.client.bucket_exists(bucket):
                    self.client.make_bucket(bucket)
                    logger.info(f"Created bucket: {bucket}")
            except S3Error as e:
                logger.error(f"Failed to create bucket {bucket}: {e}")
    
    def upload_file(
        self,
        bucket: str,
        object_name: str,
        file_path: str | Path,
        content_type: str | None = None,
    ) -> str:
        """
        Upload a file to S3.
        
        Returns:
            Object path in format bucket/object_name
        """
        try:
            self.client.fput_object(
                bucket,
                object_name,
                str(file_path),
                content_type=content_type,
            )
            logger.info(f"Uploaded {file_path} to {bucket}/{object_name}")
            return f"{bucket}/{object_name}"
        except S3Error as e:
            logger.error(f"Upload failed: {e}")
            raise
    
    def upload_data(
        self,
        bucket: str,
        object_name: str,
        data: bytes | BinaryIO,
        length: int,
        content_type: str = "application/octet-stream",
    ) -> str:
        """Upload binary data to S3."""
        try:
            if isinstance(data, bytes):
                data = BytesIO(data)
            
            self.client.put_object(
                bucket,
                object_name,
                data,
                length,
                content_type=content_type,
            )
            logger.info(f"Uploaded data to {bucket}/{object_name}")
            return f"{bucket}/{object_name}"
        except S3Error as e:
            logger.error(f"Upload failed: {e}")
            raise
    
    def download_file(
        self,
        bucket: str,
        object_name: str,
        file_path: str | Path,
    ) -> Path:
        """Download file from S3 to local path."""
        try:
            self.client.fget_object(bucket, object_name, str(file_path))
            logger.info(f"Downloaded {bucket}/{object_name} to {file_path}")
            return Path(file_path)
        except S3Error as e:
            logger.error(f"Download failed: {e}")
            raise
    
    def get_object(self, bucket: str, object_name: str) -> bytes:
        """Get object data as bytes."""
        try:
            response = self.client.get_object(bucket, object_name)
            return response.read()
        except S3Error as e:
            logger.error(f"Get object failed: {e}")
            raise
        finally:
            response.close()
            response.release_conn()
    
    def delete_object(self, bucket: str, object_name: str) -> bool:
        """Delete object from S3."""
        try:
            self.client.remove_object(bucket, object_name)
            logger.info(f"Deleted {bucket}/{object_name}")
            return True
        except S3Error as e:
            logger.error(f"Delete failed: {e}")
            return False
    
    def object_exists(self, bucket: str, object_name: str) -> bool:
        """Check if object exists."""
        try:
            self.client.stat_object(bucket, object_name)
            return True
        except S3Error:
            return False
    
    def get_object_info(self, bucket: str, object_name: str) -> dict:
        """Get object metadata."""
        try:
            stat = self.client.stat_object(bucket, object_name)
            return {
                "bucket": bucket,
                "object_name": object_name,
                "size": stat.size,
                "etag": stat.etag,
                "content_type": stat.content_type,
                "last_modified": stat.last_modified,
            }
        except S3Error as e:
            logger.error(f"Stat failed: {e}")
            raise
    
    def list_objects(
        self,
        bucket: str,
        prefix: str = "",
        recursive: bool = True,
    ) -> list[dict]:
        """List objects in bucket with optional prefix."""
        objects = []
        try:
            for obj in self.client.list_objects(bucket, prefix=prefix, recursive=recursive):
                objects.append({
                    "name": obj.object_name,
                    "size": obj.size,
                    "last_modified": obj.last_modified,
                    "etag": obj.etag,
                })
        except S3Error as e:
            logger.error(f"List objects failed: {e}")
            raise
        return objects
    
    def get_presigned_url(
        self,
        bucket: str,
        object_name: str,
        expires_hours: int = 1,
    ) -> str:
        """Get presigned URL for object download."""
        from datetime import timedelta
        
        try:
            return self.client.presigned_get_object(
                bucket,
                object_name,
                expires=timedelta(hours=expires_hours),
            )
        except S3Error as e:
            logger.error(f"Presigned URL failed: {e}")
            raise
    
    # ─────────────────────────────────────────────────────────────────────
    # Convenience methods for genomic data
    # ─────────────────────────────────────────────────────────────────────
    
    def upload_fastq(
        self,
        project_id: int,
        sample_id: int,
        file_path: Path,
        read_number: int = 1,
    ) -> str:
        """Upload FASTQ file."""
        object_name = f"{project_id}/{sample_id}/{file_path.name}"
        return self.upload_file(
            settings.minio.bucket_raw,
            object_name,
            file_path,
            content_type="application/gzip",
        )
    
    def upload_bam(
        self,
        project_id: int,
        sample_id: int,
        file_path: Path,
    ) -> str:
        """Upload BAM/CRAM file."""
        object_name = f"{project_id}/{sample_id}/{file_path.name}"
        return self.upload_file(
            settings.minio.bucket_aligned,
            object_name,
            file_path,
            content_type="application/octet-stream",
        )
    
    def upload_vcf(
        self,
        project_id: int,
        sample_id: int,
        file_path: Path,
    ) -> str:
        """Upload VCF file."""
        object_name = f"{project_id}/{sample_id}/{file_path.name}"
        return self.upload_file(
            settings.minio.bucket_variants,
            object_name,
            file_path,
            content_type="application/gzip",
        )
    
    def get_sample_files(self, project_id: int, sample_id: int) -> dict:
        """Get all files for a sample."""
        prefix = f"{project_id}/{sample_id}/"
        
        return {
            "raw": self.list_objects(settings.minio.bucket_raw, prefix),
            "aligned": self.list_objects(settings.minio.bucket_aligned, prefix),
            "variants": self.list_objects(settings.minio.bucket_variants, prefix),
            "reports": self.list_objects(settings.minio.bucket_reports, prefix),
        }


# Lazy initialization (only when MinIO is available)
_s3_service: S3StorageService | None = None


def get_s3_service() -> S3StorageService:
    """Get S3 storage service instance."""
    global _s3_service
    
    if _s3_service is None:
        _s3_service = S3StorageService()
    
    return _s3_service
