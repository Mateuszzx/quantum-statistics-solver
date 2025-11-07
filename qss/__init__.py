"""
Quantum Statistics Solver - A library for quantum statistics calculations
"""

__version__ = "0.1.0"

# Import C++ extension (if available)
try:
    from ._qss import *  # Import C++ bindings with underscore prefix
except ImportError:
    # C++ extension not built yet (e.g., during editable install without build)
    pass

# Import Python utilities
from .data_format import numpy4d_to_dataframe
from .prep import prepare_2D_lattice_data

__all__ = [
    # Python utilities
    "numpy4d_to_dataframe",
    "prepare_2D_lattice_data",
    # C++ functions will be added here when extension is imported
]
