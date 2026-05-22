"""Solve chemical potential μ from band filling *n* (standalone / Distributions view)."""

from __future__ import annotations

import numpy as np
from qss.chemical_potential import chemical_potential
from qss import solver_options as so


def solve_mu_from_n(
    *,
    e_k: np.ndarray,
    distro: str,
    n: float,
    kbt: float,
    h_a: float,
    u: float,
) -> tuple[float | None, str | None]:
    """Same Brent defaults as ``ScientificPlotter``."""
    solver_opts = (
        so.brent_defaults(
            abs_tol=1e-9,
            rel_tol=1e-6,
            expand_factor=1.4,
            initial_step=10.0,
            max_expands=64,
            max_iters=5000,
        )
    ).to_qss()
    df = chemical_potential(
        distros=[distro],
        n_vals=[float(n)],
        kBT_vals=[float(kbt)],
        H_a_vals=[float(h_a)],
        E=e_k,
        U_vals=[float(u)],
        N=int(e_k.shape[0]),
        solver_opts=solver_opts,
        x0=2.0,
        a0=-100.0,
        b0=100.0,
    )
    if df.empty:
        return None, "Solver returned no result."
    row = df.iloc[0]
    if not bool(row["converged"]):
        return None, f"Did not converge (residual = {float(row['residual']):.3e})."
    return float(row["mu"]), None
