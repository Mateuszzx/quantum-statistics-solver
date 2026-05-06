"""Backward-compatible re-export (deprecated).

Prefer importing from ``gui.paths``.
"""

from gui.paths import *  # noqa: F403

__all__ = ["ensure_qss_results_tools_on_path", "qss_results_repo_root", "quantum_statistics_solver_root"]
