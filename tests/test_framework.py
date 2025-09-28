"""
Simple test to verify test framework is working
"""
import pytest
from pathlib import Path


def test_pytest_working():
    """Verify pytest is working"""
    assert True


def test_project_root_exists(project_root):
    """Test that project root fixture works"""
    assert project_root.exists()
    assert (project_root / "src").exists()


def test_fixtures_working(sample_vcf_content, sample_fastq_content):
    """Test that fixtures are working"""
    assert "##fileformat=VCFv4.2" in sample_vcf_content
    assert "@SEQ1" in sample_fastq_content


@pytest.mark.unit
def test_unit_marker():
    """Test with unit marker"""
    assert 1 + 1 == 2


@pytest.mark.integration
def test_integration_marker():
    """Test with integration marker"""
    assert True


@pytest.mark.slow
def test_slow_marker():
    """Test with slow marker (should be skipped by default)"""
    import time
    time.sleep(0.1)
    assert True
