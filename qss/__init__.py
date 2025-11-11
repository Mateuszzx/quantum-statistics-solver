"""QSS: Quantum Statistics Solver package."""
from __future__ import annotations

# Import extension module first to avoid circular import
from . import _qss
from .chemical_potential import chemical_potential

# Then import Python wrappers
from .distributions import DistType, prepare_2D_lattice_data
from .solver_options import (
    BracketerType,
    CoolingSchedule,
    SolverType,
    brent_defaults,
    sa_defaults,
)

__all__ = [
    "BracketerType",
    "CoolingSchedule",
    "DistType",
    "SolverType",
    "_qss",
    "brent_defaults",
    "chemical_potential",
    "prepare_2D_lattice_data",
    "sa_defaults",
]
