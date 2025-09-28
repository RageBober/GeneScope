"""
Basic tests for GenoScope project structure and setup
"""
import pytest
import sys
from pathlib import Path

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


class TestProjectSetup:
    """Test project structure and configuration"""
    
    def test_project_structure(self):
        """Verify project directory structure"""
        assert (PROJECT_ROOT / "src").exists(), "src directory not found"
        assert (PROJECT_ROOT / "frontend").exists(), "frontend directory not found"
        assert (PROJECT_ROOT / "tests").exists(), "tests directory not found"
        assert (PROJECT_ROOT / "scripts").exists(), "scripts directory not found"
        assert (PROJECT_ROOT / "monitoring").exists(), "monitoring directory not found"
        assert (PROJECT_ROOT / "k8s").exists(), "k8s directory not found"
    
    def test_configuration_files(self):
        """Verify all configuration files exist"""
        config_files = [
            "docker-compose.yml",
            "Dockerfile",
            "requirements.txt",
            "pyproject.toml",
            "Makefile",
            ".github/workflows/ci-cd.yml"
        ]
        
        for config_file in config_files:
            assert (PROJECT_ROOT / config_file).exists(), f"{config_file} not found"
    
    def test_python_modules(self):
        """Verify Python module structure"""
        src_dir = PROJECT_ROOT / "src" / "genoscope"
        
        # Check main modules
        assert (src_dir / "__init__.py").exists(), "genoscope __init__.py not found"
        assert (src_dir / "api").exists(), "API module not found"
        assert (src_dir / "pipeline").exists(), "Pipeline module not found"
        assert (src_dir / "core").exists(), "Core module not found"
        
        # Check API structure
        api_dir = src_dir / "api"
        assert (api_dir / "main.py").exists(), "API main.py not found"
        assert (api_dir / "genomics_router.py").exists(), "Genomics router not found"
        assert (api_dir / "pipeline_router.py").exists(), "Pipeline router not found"
        assert (api_dir / "billing_router.py").exists(), "Billing router not found"
    
    def test_pipeline_modules(self):
        """Verify pipeline modules exist"""
        pipeline_dir = PROJECT_ROOT / "src" / "genoscope" / "pipeline"
        
        modules = [
            "qc.py",
            "alignment.py",
            "variant_calling.py",
            "main_pipeline.py"
        ]
        
        for module in modules:
            assert (pipeline_dir / module).exists(), f"Pipeline module {module} not found"
    
    def test_frontend_structure(self):
        """Verify frontend structure"""
        frontend_dir = PROJECT_ROOT / "frontend"
        
        assert (frontend_dir / "package.json").exists(), "package.json not found"
        assert (frontend_dir / "src").exists(), "Frontend src not found"
        assert (frontend_dir / "public").exists(), "Frontend public not found"
        assert (frontend_dir / "src" / "App.tsx").exists(), "App.tsx not found"
    
    def test_docker_files(self):
        """Verify Docker configuration"""
        assert (PROJECT_ROOT / "Dockerfile").exists(), "Dockerfile not found"
        assert (PROJECT_ROOT / "docker-compose.yml").exists(), "docker-compose.yml not found"
        assert (PROJECT_ROOT / "frontend" / "Dockerfile").exists(), "Frontend Dockerfile not found"
    
    def test_kubernetes_files(self):
        """Verify Kubernetes configuration"""
        k8s_dir = PROJECT_ROOT / "k8s"
        
        assert k8s_dir.exists(), "k8s directory not found"
        assert (k8s_dir / "base").exists(), "k8s base directory not found"
        assert (k8s_dir / "base" / "deployment.yaml").exists(), "deployment.yaml not found"
    
    def test_monitoring_setup(self):
        """Verify monitoring configuration"""
        monitoring_dir = PROJECT_ROOT / "monitoring"
        
        assert monitoring_dir.exists(), "monitoring directory not found"
        assert (monitoring_dir / "docker-compose.yml").exists(), "Monitoring docker-compose not found"
        assert (monitoring_dir / "prometheus.yml").exists(), "Prometheus config not found"
        assert (monitoring_dir / "alert_rules.yml").exists(), "Alert rules not found"


class TestImports:
    """Test that modules can be imported"""
    
    def test_import_basic_modules(self):
        """Test basic module imports"""
        try:
            import fastapi
            import sqlalchemy
            import pandas
            import numpy
            assert True
        except ImportError as e:
            pytest.fail(f"Failed to import basic module: {e}")
    
    def test_import_qc_module(self):
        """Test QC module import"""
        try:
            from src.genoscope.pipeline.qc import QualityController, QCMetrics
            
            # Create instance to verify it works
            metrics = QCMetrics()
            assert hasattr(metrics, 'total_reads')
            assert hasattr(metrics, 'to_dict')
            
        except ImportError as e:
            pytest.skip(f"QC module not available: {e}")
    
    def test_import_variant_module(self):
        """Test variant calling module import"""
        try:
            from src.genoscope.pipeline.variant_calling import VariantCaller, VariantStats
            
            # Create instance to verify it works
            stats = VariantStats()
            assert hasattr(stats, 'total_variants')
            assert hasattr(stats, 'snps')
            
        except ImportError as e:
            pytest.skip(f"Variant module not available: {e}")
    
    def test_import_alignment_module(self):
        """Test alignment module import"""
        try:
            from src.genoscope.pipeline.alignment import AlignmentEngine, AlignmentStats
            
            # Create instance to verify it works
            stats = AlignmentStats()
            assert hasattr(stats, 'total_reads')
            assert hasattr(stats, 'mapped_reads')
            
        except ImportError as e:
            pytest.skip(f"Alignment module not available: {e}")


class TestDependencies:
    """Test external dependencies"""
    
    def test_required_packages(self):
        """Test that all required packages are installed"""
        required_packages = [
            "fastapi",
            "uvicorn",
            "pandas",
            "numpy",
            "sqlalchemy",
            "pydantic",
            "pytest",
            "redis",
            "celery",
        ]
        
        for package in required_packages:
            try:
                __import__(package)
            except ImportError:
                pytest.fail(f"Required package not installed: {package}")
    
    def test_optional_packages(self):
        """Test optional packages (skip if not installed)"""
        optional_packages = [
            "locust",
            "cypress",
            "stripe",
            "boto3",
        ]
        
        for package in optional_packages:
            try:
                __import__(package)
                print(f"✓ {package} installed")
            except ImportError:
                print(f"✗ {package} not installed (optional)")
