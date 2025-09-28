"""
Cloud Storage Integration Module for GenoScope

Provides unified interface for cloud storage providers (S3, Azure, GCS).
"""

from __future__ import annotations

import hashlib
import logging
import os
from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timedelta
from io import BytesIO
from pathlib import Path
from typing import Any, BinaryIO, Dict, Iterator, List, Optional, Tuple, Union
from urllib.parse import urlparse

logger = logging.getLogger(__name__)


@dataclass
class StorageObject:
    """Represents a file/object in cloud storage."""
    
    key: str
    size: int
    last_modified: datetime
    etag: Optional[str] = None
    storage_class: Optional[str] = None
    metadata: Dict[str, str] = None
    
    def __post_init__(self):
        if self.metadata is None:
            self.metadata = {}


@dataclass
class UploadResult:
    """Result of an upload operation."""
    
    success: bool
    key: str
    size: int
    etag: Optional[str] = None
    error: Optional[str] = None
    duration: Optional[float] = None


class CloudStorageInterface(ABC):
    """Abstract base class for cloud storage providers."""
    
    @abstractmethod
    def upload_file(
        self,
        local_path: Path,
        remote_key: str,
        metadata: Optional[Dict[str, str]] = None
    ) -> UploadResult:
        """Upload a file to cloud storage."""
        pass
    
    @abstractmethod
    def download_file(
        self,
        remote_key: str,
        local_path: Path
    ) -> bool:
        """Download a file from cloud storage."""
        pass
    
    @abstractmethod
    def stream_download(
        self,
        remote_key: str,
        chunk_size: int = 8192
    ) -> Iterator[bytes]:
        """Stream download a file."""
        pass
    
    @abstractmethod
    def delete_file(self, remote_key: str) -> bool:
        """Delete a file from cloud storage."""
        pass
    
    @abstractmethod
    def list_objects(
        self,
        prefix: Optional[str] = None,
        max_results: Optional[int] = None
    ) -> List[StorageObject]:
        """List objects in storage."""
        pass
    
    @abstractmethod
    def get_presigned_url(
        self,
        remote_key: str,
        expiration: int = 3600,
        operation: str = "get"
    ) -> str:
        """Generate a presigned URL."""
        pass
    
    @abstractmethod
    def file_exists(self, remote_key: str) -> bool:
        """Check if a file exists."""
        pass
    
    @abstractmethod
    def get_file_info(self, remote_key: str) -> Optional[StorageObject]:
        """Get file metadata."""
        pass


class S3Storage(CloudStorageInterface):
    """Amazon S3 storage implementation."""
    
    def __init__(
        self,
        bucket_name: str,
        region: str = "us-east-1",
        access_key: Optional[str] = None,
        secret_key: Optional[str] = None,
        endpoint_url: Optional[str] = None  # For S3-compatible services
    ):
        """
        Initialize S3 storage client.
        
        Args:
            bucket_name: S3 bucket name
            region: AWS region
            access_key: AWS access key (optional, uses env if not provided)
            secret_key: AWS secret key (optional, uses env if not provided)
            endpoint_url: Custom endpoint URL for S3-compatible services
        """
        try:
            import boto3
            from botocore.exceptions import NoCredentialsError
        except ImportError:
            raise ImportError("boto3 is required for S3 storage. Install with: pip install boto3")
        
        self.bucket_name = bucket_name
        self.region = region
        
        # Create S3 client
        session_kwargs = {}
        if access_key and secret_key:
            session_kwargs["aws_access_key_id"] = access_key
            session_kwargs["aws_secret_access_key"] = secret_key
        
        client_kwargs = {"region_name": region}
        if endpoint_url:
            client_kwargs["endpoint_url"] = endpoint_url
        
        session = boto3.Session(**session_kwargs)
        self.s3_client = session.client("s3", **client_kwargs)
        self.s3_resource = session.resource("s3", **client_kwargs)
        self.bucket = self.s3_resource.Bucket(bucket_name)
        
        # Test connection
        try:
            self.s3_client.head_bucket(Bucket=bucket_name)
            logger.info(f"Successfully connected to S3 bucket: {bucket_name}")
        except Exception as e:
            logger.error(f"Failed to connect to S3 bucket {bucket_name}: {e}")
            raise
    
    def upload_file(
        self,
        local_path: Path,
        remote_key: str,
        metadata: Optional[Dict[str, str]] = None
    ) -> UploadResult:
        """Upload a file to S3."""
        import time
        from botocore.exceptions import ClientError
        
        start_time = time.time()
        
        try:
            # Prepare upload arguments
            extra_args = {}
            if metadata:
                extra_args["Metadata"] = metadata
            
            # Calculate MD5 for integrity
            with open(local_path, "rb") as f:
                file_hash = hashlib.md5(f.read()).hexdigest()
                extra_args["Metadata"] = extra_args.get("Metadata", {})
                extra_args["Metadata"]["md5"] = file_hash
            
            # Upload file
            self.s3_client.upload_file(
                str(local_path),
                self.bucket_name,
                remote_key,
                ExtraArgs=extra_args
            )
            
            # Get uploaded file info
            response = self.s3_client.head_object(
                Bucket=self.bucket_name,
                Key=remote_key
            )
            
            duration = time.time() - start_time
            
            return UploadResult(
                success=True,
                key=remote_key,
                size=response["ContentLength"],
                etag=response.get("ETag", "").strip('"'),
                duration=duration
            )
            
        except ClientError as e:
            logger.error(f"Failed to upload {local_path} to S3: {e}")
            return UploadResult(
                success=False,
                key=remote_key,
                size=0,
                error=str(e)
            )
    
    def download_file(self, remote_key: str, local_path: Path) -> bool:
        """Download a file from S3."""
        from botocore.exceptions import ClientError
        
        try:
            local_path.parent.mkdir(parents=True, exist_ok=True)
            self.s3_client.download_file(
                self.bucket_name,
                remote_key,
                str(local_path)
            )
            logger.info(f"Downloaded {remote_key} to {local_path}")
            return True
            
        except ClientError as e:
            logger.error(f"Failed to download {remote_key} from S3: {e}")
            return False
    
    def stream_download(
        self,
        remote_key: str,
        chunk_size: int = 8192
    ) -> Iterator[bytes]:
        """Stream download from S3."""
        try:
            response = self.s3_client.get_object(
                Bucket=self.bucket_name,
                Key=remote_key
            )
            
            for chunk in response["Body"].iter_chunks(chunk_size):
                yield chunk
                
        except Exception as e:
            logger.error(f"Failed to stream {remote_key}: {e}")
            raise
    
    def delete_file(self, remote_key: str) -> bool:
        """Delete a file from S3."""
        try:
            self.s3_client.delete_object(
                Bucket=self.bucket_name,
                Key=remote_key
            )
            logger.info(f"Deleted {remote_key} from S3")
            return True
            
        except Exception as e:
            logger.error(f"Failed to delete {remote_key}: {e}")
            return False
    
    def list_objects(
        self,
        prefix: Optional[str] = None,
        max_results: Optional[int] = None
    ) -> List[StorageObject]:
        """List objects in S3 bucket."""
        objects = []
        
        try:
            paginator = self.s3_client.get_paginator("list_objects_v2")
            
            page_iterator = paginator.paginate(
                Bucket=self.bucket_name,
                Prefix=prefix or "",
                PaginationConfig={"MaxItems": max_results} if max_results else {}
            )
            
            for page in page_iterator:
                if "Contents" in page:
                    for obj in page["Contents"]:
                        objects.append(StorageObject(
                            key=obj["Key"],
                            size=obj["Size"],
                            last_modified=obj["LastModified"],
                            etag=obj.get("ETag", "").strip('"'),
                            storage_class=obj.get("StorageClass")
                        ))
                        
                        if max_results and len(objects) >= max_results:
                            return objects[:max_results]
            
            return objects
            
        except Exception as e:
            logger.error(f"Failed to list objects: {e}")
            return []
    
    def get_presigned_url(
        self,
        remote_key: str,
        expiration: int = 3600,
        operation: str = "get"
    ) -> str:
        """Generate a presigned URL for S3."""
        try:
            if operation == "get":
                url = self.s3_client.generate_presigned_url(
                    "get_object",
                    Params={"Bucket": self.bucket_name, "Key": remote_key},
                    ExpiresIn=expiration
                )
            elif operation == "put":
                url = self.s3_client.generate_presigned_url(
                    "put_object",
                    Params={"Bucket": self.bucket_name, "Key": remote_key},
                    ExpiresIn=expiration
                )
            else:
                raise ValueError(f"Unsupported operation: {operation}")
            
            return url
            
        except Exception as e:
            logger.error(f"Failed to generate presigned URL: {e}")
            raise
    
    def file_exists(self, remote_key: str) -> bool:
        """Check if file exists in S3."""
        try:
            self.s3_client.head_object(
                Bucket=self.bucket_name,
                Key=remote_key
            )
            return True
        except:
            return False
    
    def get_file_info(self, remote_key: str) -> Optional[StorageObject]:
        """Get file metadata from S3."""
        try:
            response = self.s3_client.head_object(
                Bucket=self.bucket_name,
                Key=remote_key
            )
            
            return StorageObject(
                key=remote_key,
                size=response["ContentLength"],
                last_modified=response["LastModified"],
                etag=response.get("ETag", "").strip('"'),
                metadata=response.get("Metadata", {})
            )
            
        except Exception as e:
            logger.error(f"Failed to get file info for {remote_key}: {e}")
            return None


class StorageManager:
    """
    Unified storage manager for multiple providers.
    
    Manages storage across different cloud providers with
    automatic failover and load balancing.
    """
    
    def __init__(self, primary: CloudStorageInterface, secondary: Optional[CloudStorageInterface] = None):
        """
        Initialize storage manager.
        
        Args:
            primary: Primary storage provider
            secondary: Secondary storage provider for failover
        """
        self.primary = primary
        self.secondary = secondary
        self.stats = {
            "uploads": 0,
            "downloads": 0,
            "bytes_uploaded": 0,
            "bytes_downloaded": 0,
            "errors": 0
        }
    
    def upload_file(
        self,
        local_path: Path,
        remote_key: str,
        metadata: Optional[Dict[str, str]] = None,
        replicate: bool = False
    ) -> UploadResult:
        """
        Upload file with optional replication.
        
        Args:
            local_path: Local file path
            remote_key: Remote storage key
            metadata: File metadata
            replicate: Whether to replicate to secondary storage
        """
        # Upload to primary
        result = self.primary.upload_file(local_path, remote_key, metadata)
        
        if result.success:
            self.stats["uploads"] += 1
            self.stats["bytes_uploaded"] += result.size
            
            # Replicate to secondary if requested
            if replicate and self.secondary:
                secondary_result = self.secondary.upload_file(local_path, remote_key, metadata)
                if not secondary_result.success:
                    logger.warning(f"Failed to replicate {remote_key} to secondary storage")
        else:
            self.stats["errors"] += 1
            
            # Try secondary as fallback
            if self.secondary:
                logger.info("Trying secondary storage as fallback")
                result = self.secondary.upload_file(local_path, remote_key, metadata)
                if result.success:
                    self.stats["uploads"] += 1
                    self.stats["bytes_uploaded"] += result.size
        
        return result
    
    def download_file(
        self,
        remote_key: str,
        local_path: Path,
        try_secondary: bool = True
    ) -> bool:
        """
        Download file with automatic failover.
        
        Args:
            remote_key: Remote storage key
            local_path: Local file path
            try_secondary: Whether to try secondary on failure
        """
        # Try primary first
        success = self.primary.download_file(remote_key, local_path)
        
        if success:
            self.stats["downloads"] += 1
            if local_path.exists():
                self.stats["bytes_downloaded"] += local_path.stat().st_size
        else:
            self.stats["errors"] += 1
            
            # Try secondary as fallback
            if try_secondary and self.secondary:
                logger.info("Trying secondary storage for download")
                success = self.secondary.download_file(remote_key, local_path)
                if success:
                    self.stats["downloads"] += 1
                    if local_path.exists():
                        self.stats["bytes_downloaded"] += local_path.stat().st_size
        
        return success
    
    def get_stats(self) -> Dict[str, Any]:
        """Get storage manager statistics."""
        return self.stats.copy()
    
    def sync_storages(self, prefix: Optional[str] = None) -> Tuple[int, int]:
        """
        Sync files between primary and secondary storage.
        
        Args:
            prefix: Optional prefix to filter files
            
        Returns:
            Tuple of (synced_count, error_count)
        """
        if not self.secondary:
            logger.warning("No secondary storage configured for sync")
            return 0, 0
        
        synced = 0
        errors = 0
        
        # Get file lists from both storages
        primary_files = {obj.key: obj for obj in self.primary.list_objects(prefix)}
        secondary_files = {obj.key: obj for obj in self.secondary.list_objects(prefix)}
        
        # Find files missing in secondary
        missing_in_secondary = set(primary_files.keys()) - set(secondary_files.keys())
        
        for key in missing_in_secondary:
            # Download from primary and upload to secondary
            temp_path = Path(f"/tmp/{key.replace('/', '_')}")
            
            try:
                if self.primary.download_file(key, temp_path):
                    result = self.secondary.upload_file(temp_path, key)
                    if result.success:
                        synced += 1
                    else:
                        errors += 1
                else:
                    errors += 1
            except Exception as e:
                logger.error(f"Failed to sync {key}: {e}")
                errors += 1
            finally:
                # Clean up temp file
                if temp_path.exists():
                    temp_path.unlink()
        
        logger.info(f"Storage sync completed: {synced} synced, {errors} errors")
        return synced, errors
