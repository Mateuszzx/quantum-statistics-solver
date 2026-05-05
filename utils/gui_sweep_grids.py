"""Parameter sweeps (linear / log) aligned with the Jupyter notebooks."""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

SweepMode = Literal["linear", "log"]


def build_1d_sweep(
    vmin: float,
    vmax: float,
    count: int,
    mode: SweepMode,
) -> NDArray[np.float64]:
    """Return ``count`` values from ``vmin`` to ``vmax`` (inclusive along the grid)."""
    if int(count) < 2:
        return np.array([float(vmin)], dtype=np.float64)
    n = int(count)
    if mode == "linear":
        return np.linspace(float(vmin), float(vmax), n, dtype=np.float64)
    lo, hi = float(vmin), float(vmax)
    if lo <= 0 or hi <= 0:
        raise ValueError("Log sweep requires strictly positive min and max.")
    return np.logspace(np.log10(lo), np.log10(hi), n, dtype=np.float64)


def betas_for_scientific_plot_n(kBT_vals: NDArray[np.floating]) -> NDArray[np.float64]:
    """
    ``utils.scientific_plot.plot_result_n`` expects ``betas`` with ``k_BT = 1/beta`` per panel.
    """
    k = np.asarray(kBT_vals, dtype=np.float64)
    with np.errstate(divide="ignore"):
        out = 1.0 / k
    return out
