"""Parameter sweeps (linear / log) aligned with the Jupyter notebooks."""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

SweepMode = Literal["linear", "log"]

MIN_SWEEP_COUNT = 2


def parse_float_list(text: str) -> list[float]:
    """Parse comma- or whitespace-separated floats."""
    raw = [p.strip() for p in text.replace(";", ",").split(",") if p.strip()]
    if not raw:
        return []
    return [float(p) for p in raw]


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


def values_from_input(
    mode: str,
    *,
    scalar: float,
    list_text: str,
    vmin: float,
    vmax: float,
    count: int,
) -> NDArray[np.float64]:
    """
    Build a 1D parameter array from GUI input mode.

    ``mode`` is one of ``"Scalar"``, ``"List"``, ``"Linspace"``, ``"Logspace"``.
    """
    if mode == "Scalar":
        return np.array([float(scalar)], dtype=np.float64)
    if mode == "List":
        parsed = parse_float_list(list_text)
        if not parsed:
            raise ValueError("List mode needs at least one numeric value.")
        return np.array(parsed, dtype=np.float64)
    if mode == "Linspace":
        return build_1d_sweep(vmin, vmax, int(count), "linear")
    if mode == "Logspace":
        return build_1d_sweep(vmin, vmax, int(count), "log")
    raise ValueError(f"Unknown sweep input mode: {mode}")
