"""Repository layout helpers."""

from __future__ import annotations

from pathlib import Path


def quantum_statistics_solver_root() -> Path:
    """``.../quantum-statistics-solver``."""
    return Path(__file__).resolve().parents[1]
