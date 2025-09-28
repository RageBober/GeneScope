"""
Unit Tests for Variant Calling Module
"""

import pytest
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch
import sys

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.genoscope.pipeline.variant_calling import VariantCaller, VariantStats, JointGenotyping


class TestVariantStats:
    """Tests for VariantStats dataclass"""
    
    def test_variant_stats_initialization(self):
        """Test VariantStats initialization with default values"""
        stats = VariantStats()
        
        assert stats.total_variants == 0
        assert stats.snps == 0
        assert stats.indels == 0
        assert stats.insertions == 0
        assert stats.deletions == 0
        assert stats.multi_allelic == 0
        assert stats.homozygous == 0
        assert stats.heterozygous == 0
        assert stats.pass_filter == 0
        assert stats.filtered == 0
        assert stats.ti_tv_ratio == 0.0
        assert stats.mean_quality == 0.0
        assert stats.mean_depth == 0.0
    
    def test_variant_stats_custom_values(self):
        """Test VariantStats with custom values"""
        stats = VariantStats(
            total_variants=100,
            snps=80,
            indels=20,
            ti_tv_ratio=2.5
        )
        
        assert stats.total_variants == 100
        assert stats.snps == 80
        assert stats.indels == 20
        assert stats.ti_tv_ratio == 2.5


class TestVariantCaller:
    """Tests for VariantCaller class"""
    
    @pytest.fixture
    def variant_caller(self):
        """Create VariantCaller instance for testing"""
        with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as ref:
            ref.write(b">chr1\nATCGATCGATCGATCGATCG\n")
            ref_path = ref.name
        
        with tempfile.TemporaryDirectory() as tmpdir:
            caller = VariantCaller(
                reference_fasta=ref_path,
                work_dir=tmpdir,
                threads=2
            )
            yield caller
        
        # Cleanup
        Path(ref_path).unlink(missing_ok=True)
    
    def test_variant_caller_initialization(self, variant_caller):
        """Test VariantCaller initialization"""
        assert variant_caller.reference_fasta is not None
        assert variant_caller.work_dir.exists()
        assert variant_caller.threads == 2
        assert isinstance(variant_caller.callers, dict)
    
    def test_check_callers(self, variant_caller):
        """Test checking available variant callers"""
        callers = variant_caller._check_callers()
        
        assert isinstance(callers, dict)
        assert "gatk" in callers
        assert "freebayes" in callers
        assert "bcftools" in callers
        assert "strelka" in callers
        assert "deepvariant" in callers
        
        # Values should be boolean
        for caller, available in callers.items():
            assert isinstance(available, bool)
    
    def test_calculate_vcf_stats_empty_file(self, variant_caller):
        """Test calculating stats from empty VCF"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            vcf_path = vcf.name
        
        try:
            stats = variant_caller._calculate_vcf_stats(vcf_path)
            
            assert stats.total_variants == 0
            assert stats.snps == 0
            assert stats.indels == 0
            
        finally:
            Path(vcf_path).unlink()
    
    def test_calculate_vcf_stats_with_variants(self, variant_caller):
        """Test calculating stats from VCF with variants"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            # SNP
            vcf.write("chr1\t100\t.\tA\tG\t30\tPASS\tDP=20\tGT\t0/1\n")
            # Another SNP
            vcf.write("chr1\t200\t.\tC\tT\t40\tPASS\tDP=25\tGT\t1/1\n")
            # Deletion
            vcf.write("chr1\t300\t.\tATCG\tA\t35\tPASS\tDP=22\tGT\t0/1\n")
            # Insertion
            vcf.write("chr1\t400\t.\tA\tATCG\t45\tPASS\tDP=30\tGT\t1/1\n")
            vcf_path = vcf.name
        
        try:
            stats = variant_caller._calculate_vcf_stats(vcf_path)
            
            assert stats.total_variants == 4
            assert stats.snps == 2
            assert stats.indels == 2
            assert stats.insertions == 1
            assert stats.deletions == 1
            assert stats.pass_filter == 4
            assert stats.mean_quality > 0
            assert stats.mean_depth > 0
            
        finally:
            Path(vcf_path).unlink()
    
    def test_count_variants(self, variant_caller):
        """Test counting variants in VCF"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as vcf:
            vcf.write("##fileformat=VCFv4.2\n")
            vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            vcf.write("chr1\t100\t.\tA\tG\t30\tPASS\t.\n")
            vcf.write("chr1\t200\t.\tC\tT\t40\tPASS\t.\n")
            vcf.write("chr1\t300\t.\tG\tA\t35\tPASS\t.\n")
            vcf_path = vcf.name
        
        try:
            count = variant_caller._count_variants(vcf_path)
            assert count == 3
            
        finally:
            Path(vcf_path).unlink()


class TestJointGenotyping:
    """Tests for JointGenotyping class"""
    
    @pytest.fixture
    def joint_genotyper(self):
        """Create JointGenotyping instance for testing"""
        with tempfile.NamedTemporaryFile(suffix='.fa', delete=False) as ref:
            ref.write(b">chr1\nATCGATCGATCGATCGATCG\n")
            ref_path = ref.name
        
        with tempfile.TemporaryDirectory() as tmpdir:
            genotyper = JointGenotyping(
                reference_fasta=ref_path,
                work_dir=tmpdir
            )
            yield genotyper
        
        # Cleanup
        Path(ref_path).unlink(missing_ok=True)
    
    def test_joint_genotyping_initialization(self, joint_genotyper):
        """Test JointGenotyping initialization"""
        assert joint_genotyper.reference_fasta is not None
        assert joint_genotyper.work_dir.exists()
    
    @patch('subprocess.run')
    def test_create_genomicsdb_success(self, mock_run, joint_genotyper):
        """Test successful GenomicsDB creation"""
        mock_run.return_value = MagicMock(returncode=0)
        
        gvcf_files = [
            "/path/to/sample1.g.vcf",
            "/path/to/sample2.g.vcf",
            "/path/to/sample3.g.vcf"
        ]
        
        success = joint_genotyper.create_genomicsdb(
            gvcf_files=gvcf_files,
            interval="chr1:1-1000000",
            db_path="/tmp/genomicsdb"
        )
        
        assert success is True
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_create_genomicsdb_failure(self, mock_run, joint_genotyper):
        """Test GenomicsDB creation failure"""
        mock_run.side_effect = Exception("GATK not found")
        
        gvcf_files = ["/path/to/sample1.g.vcf"]
        
        success = joint_genotyper.create_genomicsdb(
            gvcf_files=gvcf_files,
            interval="chr1:1-1000000",
            db_path="/tmp/genomicsdb"
        )
        
        assert success is False
