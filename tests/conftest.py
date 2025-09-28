"""
Pytest configuration for GenoScope tests
This file configures paths and fixtures for all tests
"""
import sys
import os
from pathlib import Path
import pytest

# Get the project root directory
PROJECT_ROOT = Path(__file__).parent.parent
SRC_DIR = PROJECT_ROOT / "src"

# Add project directories to Python path
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(SRC_DIR))

# Set environment variables for testing
os.environ["TESTING"] = "true"
os.environ["DATABASE_URL"] = "sqlite:///./test.db"
os.environ["REDIS_URL"] = "redis://localhost:6379/1"

# Configure pytest
def pytest_configure(config):
    """Configure pytest with custom markers"""
    config.addinivalue_line("markers", "unit: Unit tests")
    config.addinivalue_line("markers", "integration: Integration tests")
    config.addinivalue_line("markers", "e2e: End-to-end tests")
    config.addinivalue_line("markers", "slow: Slow tests")
    config.addinivalue_line("markers", "performance: Performance tests")

# Shared fixtures
@pytest.fixture(scope="session")
def project_root():
    """Return project root directory"""
    return PROJECT_ROOT

@pytest.fixture(scope="session")
def test_data_dir():
    """Return test data directory"""
    test_data = PROJECT_ROOT / "tests" / "data"
    test_data.mkdir(exist_ok=True, parents=True)
    return test_data

@pytest.fixture
def temp_dir(tmp_path):
    """Create a temporary directory for test files"""
    return tmp_path

@pytest.fixture
def sample_vcf_content():
    """Sample VCF file content for testing"""
    return """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1
chr1\t100\t.\tA\tG\t30\tPASS\tDP=20\tGT:GQ\t0/1:99
chr1\t200\t.\tC\tT\t40\tPASS\tDP=25\tGT:GQ\t1/1:99
chr1\t300\t.\tATCG\tA\t35\tPASS\tDP=22\tGT:GQ\t0/1:99
"""

@pytest.fixture
def sample_fastq_content():
    """Sample FASTQ file content for testing"""
    return """@SEQ1
ATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@SEQ2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA
+
HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
@SEQ3
NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""

@pytest.fixture
def mock_config():
    """Mock configuration for testing"""
    return {
        "database_url": "sqlite:///./test.db",
        "redis_url": "redis://localhost:6379/1",
        "secret_key": "test-secret-key",
        "debug": True,
        "testing": True
    }

# Skip markers for conditional tests
def pytest_collection_modifyitems(config, items):
    """Add skip markers based on conditions"""
    
    # Skip integration tests if --unit-only flag is used
    if config.getoption("--unit-only", default=False):
        skip_integration = pytest.mark.skip(reason="Skipping integration tests")
        for item in items:
            if "integration" in item.keywords:
                item.add_marker(skip_integration)
    
    # Skip slow tests unless --include-slow flag is used
    if not config.getoption("--include-slow", default=False):
        skip_slow = pytest.mark.skip(reason="Skipping slow tests")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)

def pytest_addoption(parser):
    """Add custom command line options"""
    parser.addoption(
        "--unit-only",
        action="store_true",
        default=False,
        help="Run only unit tests"
    )
    parser.addoption(
        "--include-slow",
        action="store_true",
        default=False,
        help="Include slow tests"
    )
    parser.addoption(
        "--include-performance",
        action="store_true",
        default=False,
        help="Include performance tests"
    )
