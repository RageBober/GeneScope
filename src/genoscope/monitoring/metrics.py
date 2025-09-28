"""
Prometheus metrics integration for GenoScope
"""

from prometheus_client import Counter, Histogram, Gauge, generate_latest, CONTENT_TYPE_LATEST
from prometheus_client import CollectorRegistry
from fastapi import Request, Response
from fastapi.responses import PlainTextResponse
import time
import psutil
import logging

logger = logging.getLogger(__name__)

# Create metrics registry
registry = CollectorRegistry()

# Request metrics
request_count = Counter(
    'fastapi_requests_total',
    'Total number of requests',
    ['method', 'endpoint', 'status'],
    registry=registry
)

request_duration = Histogram(
    'fastapi_request_duration_seconds',
    'Request duration in seconds',
    ['method', 'endpoint'],
    registry=registry
)

# Pipeline metrics
pipeline_jobs_total = Counter(
    'pipeline_jobs_total',
    'Total number of pipeline jobs',
    ['status'],
    registry=registry
)

pipeline_jobs_duration = Histogram(
    'pipeline_job_duration_seconds',
    'Pipeline job duration in seconds',
    ['sample_type', 'status'],
    registry=registry
)

pipeline_jobs_active = Gauge(
    'pipeline_jobs_active',
    'Number of active pipeline jobs',
    registry=registry
)

# Variant metrics
variants_processed = Counter(
    'variants_processed_total',
    'Total number of variants processed',
    ['variant_type'],
    registry=registry
)

variant_quality_score = Histogram(
    'variant_quality_score',
    'Distribution of variant quality scores',
    registry=registry
)

# File metrics
files_uploaded = Counter(
    'files_uploaded_total',
    'Total number of files uploaded',
    ['file_type'],
    registry=registry
)

file_size_bytes = Histogram(
    'file_size_bytes',
    'Size of uploaded files in bytes',
    ['file_type'],
    registry=registry
)

# Database metrics
database_connections = Gauge(
    'database_connections_active',
    'Number of active database connections',
    registry=registry
)

database_query_duration = Histogram(
    'database_query_duration_seconds',
    'Database query duration',
    ['query_type'],
    registry=registry
)

# Billing metrics
stripe_payment_total = Counter(
    'stripe_payment_total',
    'Total number of Stripe payments',
    ['status'],
    registry=registry
)

stripe_payment_amount = Histogram(
    'stripe_payment_amount_dollars',
    'Payment amounts in dollars',
    registry=registry
)

subscription_count = Gauge(
    'subscription_count',
    'Number of active subscriptions',
    ['plan'],
    registry=registry
)

# User metrics
user_usage_percentage = Gauge(
    'user_usage_percentage',
    'User usage as percentage of plan limits',
    ['user_id', 'metric_type'],
    registry=registry
)

# System metrics
system_cpu_percent = Gauge(
    'system_cpu_percent',
    'System CPU usage percentage',
    registry=registry
)

system_memory_percent = Gauge(
    'system_memory_percent',
    'System memory usage percentage',
    registry=registry
)

system_disk_usage_percent = Gauge(
    'system_disk_usage_percent',
    'System disk usage percentage',
    ['mount_point'],
    registry=registry
)

# Celery metrics
celery_queue_length = Gauge(
    'celery_queue_length',
    'Number of tasks in Celery queue',
    ['queue_name'],
    registry=registry
)

celery_task_duration = Histogram(
    'celery_task_duration_seconds',
    'Celery task execution duration',
    ['task_name'],
    registry=registry
)

celery_task_total = Counter(
    'celery_task_total',
    'Total number of Celery tasks',
    ['task_name', 'status'],
    registry=registry
)


class MetricsMiddleware:
    """Middleware to collect request metrics"""
    
    def __init__(self, app):
        self.app = app
    
    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        
        start_time = time.time()
        
        async def send_wrapper(message):
            if message["type"] == "http.response.start":
                duration = time.time() - start_time
                
                # Record metrics
                path = scope["path"]
                method = scope["method"]
                status = message.get("status", 0)
                
                request_count.labels(
                    method=method,
                    endpoint=path,
                    status=status
                ).inc()
                
                request_duration.labels(
                    method=method,
                    endpoint=path
                ).observe(duration)
            
            await send(message)
        
        await self.app(scope, receive, send_wrapper)


def update_system_metrics():
    """Update system resource metrics"""
    try:
        # CPU usage
        system_cpu_percent.set(psutil.cpu_percent(interval=1))
        
        # Memory usage
        memory = psutil.virtual_memory()
        system_memory_percent.set(memory.percent)
        
        # Disk usage
        for partition in psutil.disk_partitions():
            if partition.mountpoint in ['/', '/data']:
                usage = psutil.disk_usage(partition.mountpoint)
                system_disk_usage_percent.labels(
                    mount_point=partition.mountpoint
                ).set(usage.percent)
    except Exception as e:
        logger.error(f"Failed to update system metrics: {e}")


async def metrics_endpoint(request: Request) -> Response:
    """Endpoint to expose metrics for Prometheus"""
    update_system_metrics()
    
    # Generate metrics in Prometheus format
    metrics = generate_latest(registry)
    
    return Response(
        content=metrics,
        media_type=CONTENT_TYPE_LATEST
    )


def record_pipeline_job(status: str, duration: float = None, sample_type: str = "wgs"):
    """Record pipeline job metrics"""
    pipeline_jobs_total.labels(status=status).inc()
    
    if duration:
        pipeline_jobs_duration.labels(
            sample_type=sample_type,
            status=status
        ).observe(duration)


def record_variant_metrics(variant_type: str, quality: float):
    """Record variant processing metrics"""
    variants_processed.labels(variant_type=variant_type).inc()
    variant_quality_score.observe(quality)


def record_file_upload(file_type: str, size_bytes: int):
    """Record file upload metrics"""
    files_uploaded.labels(file_type=file_type).inc()
    file_size_bytes.labels(file_type=file_type).observe(size_bytes)


def record_payment(status: str, amount: float = None):
    """Record payment metrics"""
    stripe_payment_total.labels(status=status).inc()
    
    if amount and status == "success":
        stripe_payment_amount.observe(amount)


def update_subscription_count(plan_counts: dict):
    """Update subscription count by plan"""
    for plan, count in plan_counts.items():
        subscription_count.labels(plan=plan).set(count)


def record_database_query(query_type: str, duration: float):
    """Record database query metrics"""
    database_query_duration.labels(query_type=query_type).observe(duration)


def update_celery_metrics(queue_lengths: dict, active_tasks: dict):
    """Update Celery queue metrics"""
    for queue_name, length in queue_lengths.items():
        celery_queue_length.labels(queue_name=queue_name).set(length)
    
    for task_name, count in active_tasks.items():
        celery_task_total.labels(
            task_name=task_name,
            status="active"
        ).inc(count)


def record_celery_task(task_name: str, status: str, duration: float = None):
    """Record Celery task execution"""
    celery_task_total.labels(
        task_name=task_name,
        status=status
    ).inc()
    
    if duration:
        celery_task_duration.labels(task_name=task_name).observe(duration)


def record_user_usage(user_id: str, usage_data: dict):
    """Record user usage metrics"""
    for metric_type, percentage in usage_data.items():
        user_usage_percentage.labels(
            user_id=user_id,
            metric_type=metric_type
        ).set(percentage)
