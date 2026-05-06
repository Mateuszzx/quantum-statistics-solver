"""Tight-binding dispersion grids for the GUI (1D / 2D), matching notebook conventions."""

from __future__ import annotations

from typing import Literal

import numpy as np
from numpy.typing import NDArray

from qss.distributions import prepare_2D_lattice_data

Dim = Literal["1D", "2D"]


def build_dispersion(
    dim: Dim,
    N: int,
    t: float,
    t_prim: float,
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating] | None]:
    """
    Return ``(e_k, kx, ky)``.

    * 2D: square lattice ``N x N`` on ``[-pi, pi)`` (same as notebooks).
    * 1D: chain dispersion on ``kx``, ``ky`` is ``None``.
    """
    kx = np.linspace(-np.pi, np.pi, int(N), endpoint=False)
    if dim == "2D":
        ky = np.linspace(-np.pi, np.pi, int(N), endpoint=False)
        e_k, _, _ = prepare_2D_lattice_data(kx, ky, float(t), float(t_prim))
        return e_k, kx, ky
    # 1D chain: use diagonal (N, N) energies so ``E`` stays square like the 2D mesh case.
    e1 = 2.0 * float(t) * np.cos(kx) + 2.0 * float(t_prim) * np.cos(2.0 * kx)
    e_k = np.diag(e1.astype(np.float64))
    return e_k, kx, None
