"""
GenoScope Pipeline Module
Genomic data processing pipelines and analysis workflows.
"""

from .genomics_pipeline import GenomicsPipeline

__all__ = ["GenomicsPipeline"]

# Legacy PipelineExecutor for backward compatibility
class PipelineExecutor:
    """Legacy pipeline executor - use GenomicsPipeline instead."""
    
    def __init__(self):
        self.pipeline = GenomicsPipeline()
    
    def execute(self, *args, **kwargs):
        """Execute pipeline - delegates to GenomicsPipeline."""
        return self.pipeline.run_full_pipeline(*args, **kwargs)
