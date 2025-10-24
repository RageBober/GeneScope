"""
Tests for bioinformatics_tools validators module
"""

import pytest
from pathlib import Path
from bioinformatics_tools.validators import (
    validate_file_path,
    validate_sequence,
    validate_docker_image_name,
    validate_command_args,
    ValidationError,
)


class TestValidateFilePath:
    """Tests for file path validation"""

    def test_valid_file_path(self, tmp_path):
        """Test validation of valid file path"""
        test_file = tmp_path / "test.fasta"
        test_file.write_text(">seq1\nATCG\n")

        result = validate_file_path(test_file)
        assert result.exists()
        assert result.is_file()

    def test_path_traversal_attack(self):
        """Test prevention of path traversal attacks"""
        with pytest.raises(ValidationError, match="traversal"):
            validate_file_path("../../etc/passwd")

    def test_invalid_extension(self, tmp_path):
        """Test rejection of invalid file extensions"""
        test_file = tmp_path / "test.exe"
        test_file.write_text("malicious")

        with pytest.raises(ValidationError, match="extension"):
            validate_file_path(test_file)

    def test_nonexistent_file(self):
        """Test handling of nonexistent files"""
        with pytest.raises(ValidationError, match="does not exist"):
            validate_file_path("/nonexistent/file.fasta")


class TestValidateSequence:
    """Tests for sequence validation"""

    def test_valid_nucleotide_sequence(self):
        """Test validation of valid nucleotide sequence"""
        seq = validate_sequence("ATCGATCG", seq_type="nucleotide")
        assert seq == "ATCGATCG"

    def test_invalid_nucleotide_characters(self):
        """Test rejection of invalid nucleotide characters"""
        with pytest.raises(ValidationError, match="Invalid nucleotide"):
            validate_sequence("ATCXYZ", seq_type="nucleotide")

    def test_valid_protein_sequence(self):
        """Test validation of valid protein sequence"""
        seq = validate_sequence("MKTAYIAKQR", seq_type="protein")
        assert seq == "MKTAYIAKQR"

    def test_sequence_too_short(self):
        """Test rejection of too short sequences"""
        with pytest.raises(ValidationError, match="too short"):
            validate_sequence("AT", seq_type="nucleotide", min_length=10)

    def test_sequence_too_long(self):
        """Test rejection of too long sequences"""
        with pytest.raises(ValidationError, match="too long"):
            validate_sequence("A" * 1000, seq_type="nucleotide", max_length=100)


class TestValidateDockerImageName:
    """Tests for Docker image name validation"""

    def test_valid_docker_image(self):
        """Test validation of valid Docker image names"""
        assert validate_docker_image_name("ubuntu:22.04") == "ubuntu:22.04"
        assert validate_docker_image_name("biocontainers/kraken2:latest")

    def test_image_with_registry(self):
        """Test Docker image with registry"""
        image = "quay.io/biocontainers/megahit:1.2.9"
        assert validate_docker_image_name(image) == image

    def test_shell_injection_attempt(self):
        """Test prevention of shell injection in image names"""
        with pytest.raises(ValidationError, match="dangerous"):
            validate_docker_image_name("ubuntu; rm -rf /")


class TestValidateCommandArgs:
    """Tests for command argument validation"""

    def test_valid_args(self):
        """Test validation of safe command arguments"""
        args = ["--input", "test.fasta", "--output", "result.txt"]
        result = validate_command_args(args)
        assert result == args

    def test_dangerous_characters(self):
        """Test rejection of arguments with shell metacharacters"""
        with pytest.raises(ValidationError, match="dangerous"):
            validate_command_args(["--input", "test.fasta; rm -rf /"])

    def test_empty_args(self):
        """Test handling of empty argument list"""
        result = validate_command_args([])
        assert result == []
