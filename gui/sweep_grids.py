"""Parameter sweeps (linear / log) aligned with the Jupyter notebooks."""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

SweepMode = Literal["linear", "log"]

MIN_SWEEP_COUNT = 2


def build_1d_sweep(
    vmin: float,
    vmax: float,
    count: int,
    mode: SweepMode,
) -> NDArray[np.float64]:
    """Return ``count`` values from ``vmin`` to ``vmax`` (inclusive along the grid)."""
    if int(count) < MIN_SWEEP_COUNT:
        return np.array([float(vmin)], dtype=np.float64)
    n = int(count)
    if mode == "linear":
        return np.linspace(float(vmin), float(vmax), n, dtype=np.float64)
    lo, hi = float(vmin), float(vmax)
    if lo <= 0 or hi <= 0:
        raise ValueError("Log sweep requires strictly positive min and max.")
    return np.logspace(np.log10(lo), np.log10(hi), n, dtype=np.float64)
