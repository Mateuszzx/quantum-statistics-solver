"""QSS: Quantum Statistics Solver package."""

from __future__ import annotations

# Import extension module first to avoid circular import
try:
    from . import _qss  # type: ignore[attr-defined]
except ImportError as e:
    raise ImportError("The _qss extension module is not available.") from e
from .chemical_potential import chemical_potential
from .distributions import DistType, prepare_2D_lattice_data
from .energy import energy
from .entropy import entropy
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
    "energy",
    "entropy",
    "prepare_2D_lattice_data",
    "sa_defaults",
]
