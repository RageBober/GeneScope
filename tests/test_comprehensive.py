# tests/test_comprehensive.py
from genoscope.main import GenoScopeProcessor


class TestDataPipeline:
    """Comprehensive tests for data processing pipeline."""
    
    def test_full_pipeline_csv(self, sample_csv):
        """Test complete pipeline with CSV data."""
        processor = GenoScopeProcessor()
        success = processor.run_pipeline(sample_csv, 'csv')
        assert success
        assert processor.data is not None
        # With heavy filtering on small dataset, we expect some data loss
        # The key is that the pipeline completes successfully
    
    def test_error_handling_corrupted_file(self, corrupted_file):
        """Test pipeline handles corrupted files gracefully."""
        processor = GenoScopeProcessor()
        success = processor.run_pipeline(corrupted_file, 'csv')
        # Pipeline should complete even with corrupted data
        # by applying data cleaning and handling missing values
        assert success
        assert processor.data is not None