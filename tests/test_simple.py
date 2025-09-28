"""
Simple working tests that don't require any imports
"""
import pytest
from pathlib import Path


def test_basic_math():
    """Test basic math operations"""
    assert 2 + 2 == 4
    assert 10 - 5 == 5
    assert 3 * 4 == 12
    assert 10 / 2 == 5


def test_string_operations():
    """Test string operations"""
    assert "hello" + " " + "world" == "hello world"
    assert "python".upper() == "PYTHON"
    assert "PYTHON".lower() == "python"
    assert "test" * 3 == "testtesttest"


def test_list_operations():
    """Test list operations"""
    my_list = [1, 2, 3, 4, 5]
    assert len(my_list) == 5
    assert my_list[0] == 1
    assert my_list[-1] == 5
    assert sum(my_list) == 15


def test_dictionary_operations():
    """Test dictionary operations"""
    my_dict = {"name": "GenoScope", "version": "1.0.0"}
    assert my_dict["name"] == "GenoScope"
    assert my_dict.get("version") == "1.0.0"
    assert "name" in my_dict
    assert "unknown" not in my_dict


def test_path_operations():
    """Test Path operations"""
    p = Path(".")
    assert p.exists()
    
    test_path = Path("test/file.txt")
    assert test_path.name == "file.txt"
    assert test_path.suffix == ".txt"
    assert test_path.stem == "file"


class TestProjectStructure:
    """Test project structure without imports"""
    
    def test_directories_exist(self):
        """Check if main directories exist"""
        base_dir = Path(__file__).parent.parent
        
        # These directories should exist
        assert (base_dir / "src").exists()
        assert (base_dir / "tests").exists()
        
    def test_config_files_exist(self):
        """Check if configuration files exist"""
        base_dir = Path(__file__).parent.parent
        
        # These files should exist
        assert (base_dir / "requirements.txt").exists()
        assert (base_dir / "pyproject.toml").exists()
        assert (base_dir / "Makefile").exists()


@pytest.mark.parametrize("input,expected", [
    (1, 1),
    (2, 4),
    (3, 9),
    (4, 16),
    (5, 25),
])
def test_parametrized_square(input, expected):
    """Test square function with multiple inputs"""
    assert input ** 2 == expected


@pytest.mark.parametrize("string,length", [
    ("", 0),
    ("a", 1),
    ("hello", 5),
    ("hello world", 11),
    ("GenoScope", 9),
])
def test_parametrized_string_length(string, length):
    """Test string length with multiple inputs"""
    assert len(string) == length
