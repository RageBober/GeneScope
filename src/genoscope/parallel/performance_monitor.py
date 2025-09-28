"""Real-time performance monitoring for parallel genomic processing.

This module provides comprehensive monitoring of processing performance,
memory usage, and system metrics with HTML report generation.
"""

import json
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import psutil

from ..core.logging_config import get_logger

logger = get_logger(__name__)


@dataclass
class TaskMetrics:
    """Comprehensive metrics for a processing task."""
    task_name: str
    start_time: float
    end_time: float | None = None
    peak_memory_mb: float = 0.0
    avg_cpu_percent: float = 0.0
    peak_cpu_percent: float = 0.0
    processing_speed_records_per_sec: float = 0.0
    errors_count: int = 0
    warnings_count: int = 0
    total_records: int = 0
    io_read_mb: float = 0.0
    io_write_mb: float = 0.0
    efficiency_score: float = 0.0

    # Additional metadata
    worker_count: int = 1
    chunk_count: int = 0
    failed_chunks: int = 0


class PerformanceMonitor:
    """Real-time performance monitoring with advanced metrics."""

    def __init__(self, monitoring_interval: float = 1.0):
        """Initialize performance monitor.
        
        Args:
            monitoring_interval: Seconds between metric updates
        """
        self.logger = logger
        self.monitoring_interval = monitoring_interval

        # Task tracking
        self.active_tasks: dict[str, TaskMetrics] = {}
        self.completed_tasks: dict[str, TaskMetrics] = {}

        # System monitoring
        self._monitoring_active = False
        self._monitor_thread: threading.Thread | None = None

        # Process reference
        try:
            self._process = psutil.Process()
        except Exception:
            self._process = None

        self.logger.info("Performance monitor initialized")

    def start_monitoring(self,
                        task_name: str,
                        expected_duration: float | None = None,
                        worker_count: int = 1) -> str:
        """Start monitoring a new task.
        
        Args:
            task_name: Name/description of the task
            expected_duration: Expected duration in seconds (optional)
            worker_count: Number of parallel workers
            
        Returns:
            Unique task ID for tracking
        """
        task_id = f"{task_name}_{int(time.time())}_{id(self)}"

        # Initialize task metrics
        self.active_tasks[task_id] = TaskMetrics(
            task_name=task_name,
            start_time=time.time(),
            worker_count=worker_count
        )

        self.logger.info(f"Started monitoring task: {task_name} (ID: {task_id})")
        return task_id

    def update_metrics(self,
                      task_id: str,
                      records_processed: int = 0,
                      errors: int = 0,
                      warnings: int = 0,
                      chunks_completed: int = 0,
                      chunks_failed: int = 0) -> None:
        """Update metrics for an active task.
        
        Args:
            task_id: Task identifier
            records_processed: Number of records processed
            errors: Number of errors encountered
            warnings: Number of warnings encountered
            chunks_completed: Number of chunks completed
            chunks_failed: Number of chunks that failed
        """
        if task_id not in self.active_tasks:
            self.logger.warning(f"Task {task_id} not found for metric update")
            return

        task = self.active_tasks[task_id]
        current_time = time.time()
        elapsed = current_time - task.start_time

        # Update task-specific metrics
        task.total_records = records_processed
        task.errors_count += errors
        task.warnings_count += warnings
        task.chunk_count = chunks_completed
        task.failed_chunks += chunks_failed

        # Calculate processing speed
        if elapsed > 0 and records_processed > 0:
            task.processing_speed_records_per_sec = records_processed / elapsed

        # Update system metrics
        self._update_system_metrics(task)

        self.logger.debug(f"Updated metrics for task {task.task_name}: "
                         f"{records_processed} records, {errors} errors")

    def _update_system_metrics(self, task: TaskMetrics) -> None:
        """Update system-level metrics for a task."""
        if not self._process:
            return

        try:
            # Memory metrics
            memory_info = self._process.memory_info()
            current_memory_mb = memory_info.rss / 1024 / 1024
            task.peak_memory_mb = max(task.peak_memory_mb, current_memory_mb)

            # CPU metrics
            cpu_percent = self._process.cpu_percent()
            task.peak_cpu_percent = max(task.peak_cpu_percent, cpu_percent)

            # Running average of CPU (simple exponential moving average)
            if task.avg_cpu_percent == 0:
                task.avg_cpu_percent = cpu_percent
            else:
                task.avg_cpu_percent = 0.9 * task.avg_cpu_percent + 0.1 * cpu_percent

            # I/O metrics
            try:
                io_counters = self._process.io_counters()
                task.io_read_mb = io_counters.read_bytes / 1024 / 1024
                task.io_write_mb = io_counters.write_bytes / 1024 / 1024
            except (AttributeError, psutil.AccessDenied):
                # I/O counters not available on all platforms
                pass

        except Exception as e:
            self.logger.debug(f"Could not update system metrics: {e}")

    def stop_monitoring(self, task_id: str) -> TaskMetrics:
        """Stop monitoring a task and calculate final metrics.
        
        Args:
            task_id: Task identifier
            
        Returns:
            Final task metrics
        """
        if task_id not in self.active_tasks:
            raise ValueError(f"Task {task_id} not found")

        task = self.active_tasks[task_id]
        task.end_time = time.time()

        # Calculate efficiency score
        task.efficiency_score = self._calculate_efficiency_score(task)

        # Move to completed tasks
        del self.active_tasks[task_id]
        self.completed_tasks[task_id] = task

        duration = task.end_time - task.start_time
        self.logger.info(
            f"Task {task.task_name} completed in {duration:.2f}s, "
            f"efficiency: {task.efficiency_score:.1f}/100, "
            f"speed: {task.processing_speed_records_per_sec:.0f} records/sec"
        )

        return task

    def _calculate_efficiency_score(self, task: TaskMetrics) -> float:
        """Calculate efficiency score (0-100) based on multiple factors."""
        if not task.end_time:
            return 0.0

        duration = task.end_time - task.start_time
        if duration <= 0:
            return 0.0

        # Speed factor (0-100): higher is better
        speed_factor = min(100, (task.processing_speed_records_per_sec / 1000) * 100)

        # Resource utilization factor (0-100)
        cpu_factor = min(100, task.avg_cpu_percent)

        # Error rate factor (0-100): fewer errors is better
        if task.total_records > 0:
            error_rate = (task.errors_count + task.failed_chunks) / max(1, task.total_records)
        else:
            error_rate = 1.0 if task.errors_count > 0 else 0.0

        error_factor = max(0, 100 - (error_rate * 200))

        # Parallelization effectiveness (0-100)
        if task.worker_count > 1:
            parallel_effectiveness = min(100, (task.processing_speed_records_per_sec / (task.worker_count * 50)) * 100)
        else:
            parallel_effectiveness = 100

        # Weighted average
        efficiency = (
            speed_factor * 0.30 +
            cpu_factor * 0.25 +
            error_factor * 0.25 +
            parallel_effectiveness * 0.20
        )

        return min(100, max(0, efficiency))

    def get_task_stats(self, task_id: str) -> dict[str, Any]:
        """Get comprehensive statistics for a task."""
        task = self.completed_tasks.get(task_id) or self.active_tasks.get(task_id)
        if not task:
            return {"error": f"Task {task_id} not found"}

        current_time = time.time()
        duration = (task.end_time or current_time) - task.start_time
        is_running = task.end_time is None

        return {
            "task_name": task.task_name,
            "status": "running" if is_running else "completed",
            "duration_seconds": duration,
            "start_time": datetime.fromtimestamp(task.start_time).isoformat(),
            "end_time": datetime.fromtimestamp(task.end_time).isoformat() if task.end_time else None,

            # Performance metrics
            "total_records": task.total_records,
            "processing_speed_records_per_sec": task.processing_speed_records_per_sec,
            "efficiency_score": task.efficiency_score,

            # Resource usage
            "peak_memory_mb": task.peak_memory_mb,
            "avg_cpu_percent": task.avg_cpu_percent,
            "peak_cpu_percent": task.peak_cpu_percent,
            "io_read_mb": task.io_read_mb,
            "io_write_mb": task.io_write_mb,

            # Error tracking
            "errors_count": task.errors_count,
            "warnings_count": task.warnings_count,

            # Parallel processing
            "worker_count": task.worker_count,
            "chunk_count": task.chunk_count,
            "failed_chunks": task.failed_chunks,
            "chunk_success_rate": ((task.chunk_count - task.failed_chunks) / max(1, task.chunk_count)) * 100
        }

    def get_all_stats(self) -> dict[str, dict[str, Any]]:
        """Get statistics for all tasks (active and completed)."""
        all_stats = {}

        # Add completed tasks
        for task_id in self.completed_tasks:
            all_stats[task_id] = self.get_task_stats(task_id)

        # Add active tasks
        for task_id in self.active_tasks:
            all_stats[task_id] = self.get_task_stats(task_id)

        return all_stats

    def export_report(self, output_path: Path, format: str = "html") -> Path:
        """Export comprehensive performance report."""
        all_stats = self.get_all_stats()

        if format.lower() == "json":
            with open(output_path, "w") as f:
                json.dump(all_stats, f, indent=2)

        elif format.lower() == "html":
            html_content = self._generate_html_report(all_stats)
            with open(output_path, "w", encoding="utf-8") as f:
                f.write(html_content)

        else:
            raise ValueError(f"Unsupported format: {format}")

        self.logger.info(f"Performance report exported to {output_path}")
        return output_path

    def _generate_html_report(self, stats: dict[str, dict[str, Any]]) -> str:
        """Generate comprehensive HTML performance report."""

        # Calculate summary statistics
        total_tasks = len(stats)
        completed_tasks = len([s for s in stats.values() if s["status"] == "completed"])
        avg_efficiency = sum(s.get("efficiency_score", 0) for s in stats.values()) / max(1, total_tasks)
        total_records = sum(s.get("total_records", 0) for s in stats.values())
        total_duration = sum(s.get("duration_seconds", 0) for s in stats.values())

        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>🧬 BioForge Performance Report</title>
    <meta charset="utf-8">
    <style>
        body {{ 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            margin: 0; padding: 20px; background: #f8f9fa; line-height: 1.6;
        }}
        .header {{ 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white; padding: 30px; border-radius: 10px; margin-bottom: 30px;
            text-align: center;
        }}
        .summary {{ 
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px; margin-bottom: 30px;
        }}
        .summary-card {{ 
            background: white; padding: 20px; border-radius: 8px; 
            box-shadow: 0 2px 4px rgba(0,0,0,0.1); text-align: center;
        }}
        .summary-value {{ font-size: 2em; font-weight: bold; color: #667eea; }}
        .task {{ 
            background: white; border: 1px solid #e9ecef; margin: 15px 0; 
            padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .task-header {{ 
            display: flex; justify-content: space-between; align-items: center;
            margin-bottom: 15px; border-bottom: 1px solid #e9ecef; padding-bottom: 10px;
        }}
        .task-name {{ font-size: 1.3em; font-weight: bold; }}
        .status {{ 
            padding: 4px 12px; border-radius: 20px; font-size: 0.9em; font-weight: bold;
        }}
        .status-completed {{ background: #d4edda; color: #155724; }}
        .status-running {{ background: #fff3cd; color: #856404; }}
        .metrics {{ 
            display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px; margin-top: 15px;
        }}
        .metric {{ text-align: center; }}
        .metric-value {{ 
            display: block; font-size: 1.4em; font-weight: bold; color: #495057;
        }}
        .metric-label {{ 
            display: block; font-size: 0.9em; color: #6c757d; margin-top: 5px;
        }}
        .efficiency-high {{ color: #28a745; }}
        .efficiency-medium {{ color: #ffc107; }}
        .efficiency-low {{ color: #dc3545; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 BioForge Performance Report</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
    
    <div class="summary">
        <div class="summary-card">
            <div class="summary-value">{total_tasks}</div>
            <div>Total Tasks</div>
        </div>
        <div class="summary-card">
            <div class="summary-value">{completed_tasks}</div>
            <div>Completed</div>
        </div>
        <div class="summary-card">
            <div class="summary-value">{avg_efficiency:.1f}</div>
            <div>Avg Efficiency</div>
        </div>
        <div class="summary-card">
            <div class="summary-value">{total_records:,}</div>
            <div>Records Processed</div>
        </div>
        <div class="summary-card">
            <div class="summary-value">{total_duration:.1f}s</div>
            <div>Total Time</div>
        </div>
    </div>
    
    <div class="tasks">"""

        # Add individual task details
        for task_id, task_stats in stats.items():
            efficiency = task_stats.get("efficiency_score", 0)
            efficiency_class = (
                "efficiency-high" if efficiency >= 80 else
                "efficiency-medium" if efficiency >= 50 else
                "efficiency-low"
            )

            status_class = f"status-{task_stats.get('status', 'unknown')}"

            html += f"""
        <div class="task">
            <div class="task-header">
                <div class="task-name">{task_stats.get('task_name', 'Unknown Task')}</div>
                <div class="status {status_class}">{task_stats.get('status', 'unknown').upper()}</div>
            </div>
            <div class="metrics">
                <div class="metric">
                    <span class="metric-value">{task_stats.get('duration_seconds', 0):.1f}s</span>
                    <span class="metric-label">⏱️ Duration</span>
                </div>
                <div class="metric">
                    <span class="metric-value">{task_stats.get('total_records', 0):,}</span>
                    <span class="metric-label">📊 Records</span>
                </div>
                <div class="metric">
                    <span class="metric-value">{task_stats.get('processing_speed_records_per_sec', 0):.0f}</span>
                    <span class="metric-label">⚡ Speed (rec/sec)</span>
                </div>
                <div class="metric">
                    <span class="metric-value {efficiency_class}">{efficiency:.1f}%</span>
                    <span class="metric-label">📈 Efficiency</span>
                </div>
                <div class="metric">
                    <span class="metric-value">{task_stats.get('peak_memory_mb', 0):.1f}MB</span>
                    <span class="metric-label">💾 Peak Memory</span>
                </div>
                <div class="metric">
                    <span class="metric-value">{task_stats.get('avg_cpu_percent', 0):.1f}%</span>
                    <span class="metric-label">🖥️ Avg CPU</span>
                </div>
                <div class="metric">
                    <span class="metric-value">{task_stats.get('worker_count', 1)}</span>
                    <span class="metric-label">👥 Workers</span>
                </div>
                <div class="metric">
                    <span class="metric-value">{task_stats.get('errors_count', 0)}</span>
                    <span class="metric-label">❌ Errors</span>
                </div>
            </div>
        </div>"""

        html += """
    </div>
</body>
</html>"""

        return html
