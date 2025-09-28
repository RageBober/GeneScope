"""
GenoScope - Advanced toolkit for genomic data analysis and bioinformatics.
"""

__version__ = "0.2.0"
__author__ = "RageBober"

try:
    from .main import GenoScopeProcessor
    from .main import main
    __all__ = ["GenoScopeProcessor", "__version__", "main"]
except ImportError:
    # Fallback если модули недоступны
    __all__ = ["__version__"]
