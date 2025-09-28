"""
Unit Tests for Pipeline QC Module

Тестируем:
1. Парсинг FASTQ файлов
2. Расчет метрик качества
3. Проверку пороговых значений
4. Обработку ошибок
"""

import pytest
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
import tempfile
import gzip

from src.genoscope.pipeline.qc import QualityController, QCMetrics


class TestQCMetrics:
    """Тесты для класса QCMetrics"""
    
    def test_qc_metrics_initialization(self):
        """Тест: Инициализация метрик с дефолтными значениями"""
        metrics = QCMetrics()
        
        assert metrics.total_reads == 0
        assert metrics.total_bases == 0
        assert metrics.q20_bases == 0
        assert metrics.q30_bases == 0
        assert metrics.gc_content == 0.0
        assert metrics.mean_quality == 0.0
    
    def test_qc_metrics_to_dict(self):
        """Тест: Конвертация метрик в словарь"""
        metrics = QCMetrics(
            total_reads=1000000,
            total_bases=150000000,
            q20_bases=140000000,
            q30_bases=130000000,
            gc_content=45.5
        )
        
        result = metrics.to_dict()
        
        assert result["total_reads"] == 1000000
        assert result["total_bases"] == 150000000
        assert result["q20_percentage"] == pytest.approx(93.33, 0.01)
        assert result["q30_percentage"] == pytest.approx(86.67, 0.01)
        assert result["gc_content"] == 45.5
    
    def test_qc_metrics_pass_check_success(self):
        """Тест: Проверка успешного прохождения QC"""
        metrics = QCMetrics(
            total_reads=2000000,
            total_bases=300000000,
            q30_bases=260000000,  # 86.67% Q30
            n_content=2.0,
            adapter_content=5.0,
            duplication_rate=20.0
        )
        
        passed, failures = metrics.is_pass()
        
        assert passed is True
        assert len(failures) == 0
    
    def test_qc_metrics_pass_check_failure(self):
        """Тест: Проверка провала QC по нескольким критериям"""
        metrics = QCMetrics(
            total_reads=500000,  # Too few reads
            total_bases=75000000,
            q30_bases=50000000,  # 66.67% Q30 - too low
            n_content=10.0,  # Too high
            adapter_content=15.0,  # Too high
            duplication_rate=60.0  # Too high
        )
        
        passed, failures = metrics.is_pass()
        
        assert passed is False
        assert len(failures) == 5  # All criteria failed
        assert any("Low read count" in f for f in failures)
        assert any("Low Q30" in f for f in failures)
        assert any("High N content" in f for f in failures)
        assert any("High adapter content" in f for f in failures)
        assert any("High duplication" in f for f in failures)
    
    def test_qc_metrics_custom_thresholds(self):
        """Тест: Использование кастомных порогов"""
        metrics = QCMetrics(
            total_reads=1500000,
            total_bases=225000000,
            q30_bases=200000000  # 88.89% Q30
        )
        
        custom_thresholds = {
            "min_reads": 2000000,  # Higher threshold
            "min_q30": 90,  # Higher Q30 requirement
            "max_n_content": 5,
            "max_adapter": 10,
            "max_duplication": 50
        }
        
        passed, failures = metrics.is_pass(custom_thresholds)
        
        assert passed is False
        assert len(failures) == 2
        assert any("Low read count" in f for f in failures)
        assert any("Low Q30" in f for f in failures)


class TestQualityController:
    """Тесты для класса QualityController"""
    
    @pytest.fixture
    def qc_controller(self):
        """Фикстура для создания экземпляра QualityController"""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield QualityController(work_dir=tmpdir, threads=2)
    
    def test_quality_controller_initialization(self, qc_controller):
        """Тест: Инициализация QualityController"""
        assert qc_controller.work_dir.exists()
        assert qc_controller.threads == 2
        assert isinstance(qc_controller.tools, dict)
    
    @patch('shutil.which')
    def test_check_tools_all_available(self, mock_which):
        """Тест: Проверка доступности всех инструментов"""
        mock_which.return_value = '/usr/bin/tool'
        
        with tempfile.TemporaryDirectory() as tmpdir:
            qc = QualityController(work_dir=tmpdir)
            
        assert qc.tools["fastqc"] is True
        assert qc.tools["fastp"] is True
        assert qc.tools["trimmomatic"] is True
        assert qc.tools["cutadapt"] is True
        assert qc.tools["multiqc"] is True
    
    @patch('shutil.which')
    def test_check_tools_none_available(self, mock_which):
        """Тест: Ни один инструмент не доступен"""
        mock_which.return_value = None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            qc = QualityController(work_dir=tmpdir)
            
        assert qc.tools["fastqc"] is False
        assert qc.tools["fastp"] is False
        assert all(not available for available in qc.tools.values())
    
    def test_basic_qc_with_valid_fastq(self, qc_controller):
        """Тест: Базовый QC анализ валидного FASTQ файла"""
        # Создаем тестовый FASTQ файл
        fastq_content = """@SEQ1
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fastq', delete=False) as f:
            f.write(fastq_content)
            fastq_file = f.name
        
        try:
            results = qc_controller._basic_qc([fastq_file])
            
            assert "files" in results
            assert fastq_file in results["files"]
            
            metrics = results["files"][fastq_file]
            assert metrics["total_reads"] == 3
            assert metrics["total_bases"] == 96  # 32 * 3
            assert metrics["gc_content"] > 0
            assert metrics["n_content"] > 0  # We have one sequence of Ns
            assert metrics["mean_quality"] > 0
            
        finally:
            Path(fastq_file).unlink()
    
    def test_basic_qc_with_gzipped_fastq(self, qc_controller):
        """Тест: QC анализ сжатого FASTQ файла"""
        fastq_content = b"""@SEQ1
ATCGATCGATCGATCGATCGATCGATCGATCG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
"""
        
        with tempfile.NamedTemporaryFile(suffix='.fastq.gz', delete=False) as f:
            with gzip.open(f.name, 'wb') as gz:
                gz.write(fastq_content)
            fastq_file = f.name
        
        try:
            results = qc_controller._basic_qc([fastq_file])
            
            assert "files" in results
            assert fastq_file in results["files"]
            
            metrics = results["files"][fastq_file]
            assert metrics["total_reads"] == 1
            assert metrics["total_bases"] == 32
            
        finally:
            Path(fastq_file).unlink()
    
    def test_basic_qc_with_invalid_file(self, qc_controller):
        """Тест: Обработка ошибки при невалидном файле"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("This is not a FASTQ file")
            invalid_file = f.name
        
        try:
            results = qc_controller._basic_qc([invalid_file])
            
            assert "files" in results
            assert invalid_file in results["files"]
            assert "error" in results["files"][invalid_file]
            
        finally:
            Path(invalid_file).unlink()
    
    @patch('subprocess.run')
    def test_run_fastqc_success(self, mock_run, qc_controller):
        """Тест: Успешный запуск FastQC"""
        mock_run.return_value = MagicMock(returncode=0, stderr="")
        qc_controller.tools["fastqc"] = True
        
        with tempfile.NamedTemporaryFile(suffix='.fastq') as f:
            results = qc_controller.run_fastqc([f.name], str(qc_controller.work_dir))
        
        assert results["tool"] == "FastQC"
        assert "timestamp" in results
        assert "files" in results
        mock_run.assert_called_once()
    
    @patch('subprocess.run')
    def test_trim_with_fastp_success(self, mock_run, qc_controller):
        """Тест: Успешная обрезка адаптеров с fastp"""
        mock_run.return_value = MagicMock(returncode=0)
        qc_controller.tools["fastp"] = True
        
        # Создаем mock JSON отчет fastp
        fastp_report = {
            "summary": {
                "before_filtering": {
                    "total_reads": 1000000,
                    "total_bases": 150000000,
                    "q30_rate": 0.85
                },
                "after_filtering": {
                    "total_reads": 950000,
                    "total_bases": 140000000,
                    "q30_rate": 0.92
                }
            },
            "filtering_result": {
                "passed_filter_reads": 950000,
                "low_quality_reads": 30000,
                "too_many_N_reads": 20000
            }
        }
        
        with tempfile.NamedTemporaryFile(suffix='.fastq') as f1, \
             tempfile.NamedTemporaryFile(suffix='.fastq') as f2:
            
            with patch('builtins.open', create=True) as mock_open:
                mock_open.return_value.__enter__.return_value.read.return_value = str(fastp_report)
                with patch('json.load', return_value=fastp_report):
                    
                    out_r1, out_r2, metrics = qc_controller.trim_adapters(
                        f1.name, f2.name, str(qc_controller.work_dir)
                    )
            
            assert out_r1.endswith("_trimmed_R1.fastq.gz")
            assert out_r2.endswith("_trimmed_R2.fastq.gz")
            assert metrics["tool"] == "fastp"
            assert metrics["status"] == "success"
            assert metrics["before_filtering"]["total_reads"] == 1000000
            assert metrics["after_filtering"]["total_reads"] == 950000
    
    @patch('subprocess.run')
    def test_remove_duplicates_with_picard(self, mock_run, qc_controller):
        """Тест: Удаление дупликатов с Picard"""
        mock_run.return_value = MagicMock(returncode=0)
        
        with patch('shutil.which', return_value='/usr/bin/picard'):
            with tempfile.NamedTemporaryFile(suffix='.bam') as bam_file:
                # Создаем mock metrics файл
                metrics_content = """LIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\tSECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\tUNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\tREAD_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\tESTIMATED_LIBRARY_SIZE
Unknown\t0\t5000000\t0\t0\t0\t500000\t10000\t0.1\t50000000"""
                
                with patch('builtins.open', create=True) as mock_open:
                    mock_open.return_value.__enter__.return_value.__iter__.return_value = metrics_content.split('\n')
                    
                    output_bam, metrics = qc_controller.remove_duplicates(bam_file.name)
                
                assert output_bam.endswith('.dedup.bam')
                assert metrics["tool"] == "picard"
                assert metrics["status"] == "success"
                assert metrics["duplication_rate"] == 0.1
