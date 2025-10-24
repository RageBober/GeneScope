"""
Tests for BioForge integration pipelines
"""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

from bioinformatics_tools.bioforge_integration import (
    MetagenomicsPipeline,
    ProteinAnalysisPipeline,
    MutationEffectPipeline,
    MetaGraphDataSource,
    PipelineResult,
)


class TestPipelineResult:
    """Tests for PipelineResult class"""

    def test_pipeline_result_creation(self):
        """Test creating a pipeline result"""
        result = PipelineResult(
            pipeline_name="TestPipeline",
            status="success"
        )

        assert result.pipeline_name == "TestPipeline"
        assert result.status == "success"
        assert result.steps_completed == []
        assert result.steps_failed == []

    def test_pipeline_result_to_json(self, tmp_path):
        """Test saving pipeline result to JSON"""
        result = PipelineResult(
            pipeline_name="TestPipeline",
            status="success",
            steps_completed=["step1", "step2"],
            results={"key": "value"}
        )

        output_file = tmp_path / "result.json"
        result.to_json(output_file)

        assert output_file.exists()

        # Verify content
        import json
        with open(output_file) as f:
            data = json.load(f)

        assert data["pipeline_name"] == "TestPipeline"
        assert data["status"] == "success"
        assert "step1" in data["steps_completed"]


class TestMetagenomicsPipeline:
    """Tests for MetagenomicsPipeline"""

    @pytest.fixture
    def mock_databases(self, tmp_path):
        """Create mock database directories"""
        kraken2_db = tmp_path / "kraken2_db"
        gtdbtk_db = tmp_path / "gtdbtk_db"

        kraken2_db.mkdir()
        gtdbtk_db.mkdir()

        # Create dummy database files
        (kraken2_db / "hash.k2d").touch()
        (gtdbtk_db / "metadata").mkdir()

        return kraken2_db, gtdbtk_db

    def test_pipeline_initialization(self, mock_databases):
        """Test initializing metagenomics pipeline"""
        kraken2_db, gtdbtk_db = mock_databases

        pipeline = MetagenomicsPipeline(
            kraken2_db=kraken2_db,
            gtdbtk_db=gtdbtk_db,
            assembler="megahit",
            use_docker=False  # Disable Docker for testing
        )

        assert pipeline.kraken2_db == kraken2_db
        assert pipeline.gtdbtk_db == gtdbtk_db
        assert pipeline.assembler == "megahit"

    def test_invalid_assembler(self, mock_databases):
        """Test rejection of invalid assembler"""
        kraken2_db, _ = mock_databases

        with pytest.raises(ValueError, match="Unknown assembler"):
            MetagenomicsPipeline(
                kraken2_db=kraken2_db,
                assembler="invalid_assembler"
            )

    @patch('bioinformatics_tools.bioforge_integration.Kraken2Tool')
    @patch('bioinformatics_tools.bioforge_integration.MEGAHITTool')
    def test_pipeline_run_success(self, mock_megahit, mock_kraken2, mock_databases, tmp_path):
        """Test successful pipeline execution"""
        kraken2_db, _ = mock_databases

        # Mock tool results
        mock_kraken2_instance = MagicMock()
        mock_kraken2_instance.classify.return_value = {
            "results": {
                "classified_percentage": 85.5,
                "total_reads": 10000,
                "taxa_found": [
                    {"rank": "S", "name": "E. coli", "percentage": 45.0},
                    {"rank": "S", "name": "S. aureus", "percentage": 30.0},
                ]
            }
        }
        mock_kraken2.return_value = mock_kraken2_instance

        mock_megahit_instance = MagicMock()
        mock_megahit_instance.assemble.return_value = {
            "results": {
                "num_contigs": 1500,
                "n50": 2500,
                "total_length": 5000000,
                "longest_contig": 15000
            }
        }
        mock_megahit.return_value = mock_megahit_instance

        # Create test input
        test_fastq = tmp_path / "test.fastq"
        test_fastq.write_text(">read1\nATCG\n")

        # Run pipeline
        pipeline = MetagenomicsPipeline(
            kraken2_db=kraken2_db,
            assembler="megahit",
            use_docker=False
        )

        result = pipeline.run(
            fastq_r1=test_fastq,
            output_dir=tmp_path,
            skip_gtdbtk=True  # Skip GTDB-Tk for simplicity
        )

        # Assertions
        assert result.status in ["success", "partial"]
        assert "kraken2_classification" in result.steps_completed
        assert "megahit_assembly" in result.steps_completed
        assert result.results["kraken2"]["classified_percentage"] == 85.5


class TestProteinAnalysisPipeline:
    """Tests for ProteinAnalysisPipeline"""

    @patch('bioinformatics_tools.bioforge_integration.XtriMoPGLMTool')
    def test_protein_pipeline_run(self, mock_xtrimopglm, tmp_path):
        """Test protein analysis pipeline"""
        # Mock XtriMoPGLM results
        mock_tool = MagicMock()
        mock_tool.predict.return_value = {
            "results": {
                "predictions": [
                    {"protein_id": "P1", "function": "kinase"},
                    {"protein_id": "P2", "function": "transporter"},
                ]
            }
        }
        mock_xtrimopglm.return_value = mock_tool

        # Create test input
        test_fasta = tmp_path / "proteins.fasta"
        test_fasta.write_text(">P1\nMKTAYIAKQR\n>P2\nMLKGVFTR\n")

        # Run pipeline
        pipeline = ProteinAnalysisPipeline(
            device="cpu",
            use_docker=False
        )

        result = pipeline.run(
            protein_fasta=test_fasta,
            output_dir=tmp_path
        )

        # Assertions
        assert result.status == "success"
        assert "xtrimopglm_prediction" in result.steps_completed
        assert len(result.results["predictions"]) == 2


class TestMutationEffectPipeline:
    """Tests for MutationEffectPipeline"""

    @patch('bioinformatics_tools.bioforge_integration.EnformerTool')
    def test_mutation_pipeline_run(self, mock_enformer, tmp_path):
        """Test mutation effect prediction pipeline"""
        # Mock Enformer results
        mock_tool = MagicMock()
        mock_tool.predict.return_value = {
            "results": {
                "predictions": [
                    {"variant_id": "V1", "effect_score": 0.85},
                    {"variant_id": "V2", "effect_score": 0.42},
                ]
            }
        }
        mock_enformer.return_value = mock_tool

        # Create test input
        test_fasta = tmp_path / "mutations.fasta"
        test_fasta.write_text(">V1\nATCGATCG\n>V2\nGCTAGCTA\n")

        # Run pipeline
        pipeline = MutationEffectPipeline(
            model_name="enformer",
            device="cpu",
            use_docker=False
        )

        result = pipeline.run(
            sequence_fasta=test_fasta,
            output_dir=tmp_path
        )

        # Assertions
        assert result.status == "success"
        assert "enformer_prediction" in result.steps_completed
        assert len(result.results["predictions"]) == 2


class TestMetaGraphDataSource:
    """Tests for MetaGraphDataSource"""

    @patch('bioinformatics_tools.bioforge_integration.MetaGraphTool')
    def test_metagraph_search(self, mock_metagraph, tmp_path):
        """Test MetaGraph sequence search"""
        # Create mock index directory
        index_dir = tmp_path / "metagraph_index"
        index_dir.mkdir()

        # Mock MetaGraph results
        mock_tool = MagicMock()
        mock_tool.search.return_value = {
            "results": {
                "num_matches": 42,
                "matches": [
                    {"query_id": "Q1", "matched_labels": ["SRR123", "SRR456"]},
                    {"query_id": "Q2", "matched_labels": ["SRR789"]},
                ]
            }
        }
        mock_metagraph.return_value = mock_tool

        # Create test query
        test_query = tmp_path / "query.fasta"
        test_query.write_text(">Q1\nATCG\n>Q2\nGCTA\n")

        # Run search
        data_source = MetaGraphDataSource(
            index_dir=index_dir,
            use_docker=False
        )

        result = data_source.search_sequences(
            query_fasta=test_query,
            discovery_fraction=0.8
        )

        # Assertions
        assert result.status == "success"
        assert result.results["num_matches"] == 42
        assert "metagraph_search" in result.steps_completed

    @patch('bioinformatics_tools.bioforge_integration.MetaGraphTool')
    def test_metagraph_build_index(self, mock_metagraph, tmp_path):
        """Test MetaGraph index building"""
        index_dir = tmp_path / "metagraph_index"
        index_dir.mkdir()

        # Mock build results
        mock_tool = MagicMock()
        mock_tool.build_index.return_value = {
            "results": {
                "status": "success",
                "index_size_gb": 5.2
            }
        }
        mock_metagraph.return_value = mock_tool

        # Create test sequences
        test_fasta = tmp_path / "sequences.fasta"
        test_fasta.write_text(">seq1\nATCG\n")

        # Build index
        data_source = MetaGraphDataSource(
            index_dir=index_dir,
            use_docker=False
        )

        result = data_source.build_index(
            input_fasta=test_fasta,
            kmer_length=31
        )

        # Assertions
        assert result.status == "success"
        assert "metagraph_build" in result.steps_completed
