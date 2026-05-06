"""Repository layout helpers (optional ``tools/scientific_plot_v2`` from the parent ``qss-results`` repo)."""

from __future__ import annotations

import sys
from pathlib import Path


def quantum_statistics_solver_root() -> Path:
    """``.../quantum-statistics-solver``."""
    return Path(__file__).resolve().parents[1]


def qss_results_repo_root() -> Path | None:
    """
    Parent of the solver checkout when embedded as ``qss-results/quantum-statistics-solver``.

    Returns ``None`` if that layout is not detected.
    """
    cand = quantum_statistics_solver_root().parent
    if (cand / "tools" / "scientific_plot_v2.py").is_file():
        return cand
    return None


def ensure_qss_results_tools_on_path() -> Path | None:
    """Insert ``qss-results`` root on ``sys.path`` so ``import tools.scientific_plot_v2`` works."""
    root = qss_results_repo_root()
    if root is None:
        return None
    s = str(root)
    if s not in sys.path:
        sys.path.insert(0, s)
    return root

