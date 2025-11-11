"""High-level wrappers for chemical potential calculations."""
from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

from . import _qss
from .distributions import DistType as PyDistType


def _coerce_native_distro(d: _qss.DistType | PyDistType | str) -> _qss.DistType:
    """Convert distribution type to native _qss.DistType."""
    if isinstance(d, _qss.DistType):
        return d
    if isinstance(d, PyDistType):
        return getattr(_qss.DistType, d.name)
    if isinstance(d, str):
        return getattr(_qss.DistType, d.upper())
    raise TypeError(f"Unsupported dist type: {type(d)}")


def _dist_label(d: _qss.DistType | PyDistType | str) -> str:
    """Get string label for distribution type."""
    if isinstance(d, (_qss.DistType, PyDistType)):
        return d.name
    return str(d).upper()


def chemical_potential( # noqa: PLR0913
    distros: Sequence[PyDistType | _qss.DistType | str],
    n_vals: Sequence[float],
    beta_vals: Sequence[float],
    H_a_vals: Sequence[float],
    E: np.ndarray,
    N: int | None = None,
    solver_opts: _qss.CombinedSolverOptions | None = None,
    g: float = 2.0,
    muB_eV_per_T: float = 5.788381806e-5,
    a0: float = -10.0,
    b0: float = 10.0,
    x0: float = 0.0,
) -> pd.DataFrame:
    """
    Compute chemical potential for multiple parameter combinations.

    Parameters
    ----------
    distros : distribution type or sequence of distribution types
        FD, FLB, NFLB, or ND (as string, enum, or native type)
    n_vals : sequence of float
        Band filling values
    beta_vals : sequence of float
        Inverse temperature values (1/k_B T)
    H_a_vals : sequence of float
        Applied magnetic field values (Tesla)
    E : ndarray
        Energy matrix (square 2D array)
    N : int, optional
        Effective system size (defaults to E.shape[0])
    solver_opts : CombinedSolverOptions, optional
        Solver configuration
    g : float
        Landé g-factor (default: 2.0)
    muB_eV_per_T : float
        Bohr magneton in eV/T (default: 5.788381806e-5)
    a0, b0, x0 : float
        Initial bracketing and guess parameters

    Returns
    -------
    pd.DataFrame
        MultiIndex DataFrame with columns: mu, iterations, converged, residual
        Index levels: distro, n, beta, H_a
    """
    if isinstance(distros, (str, _qss.DistType, PyDistType)):
        distro_iter = [distros]
    else:
        distro_iter = list(distros)

    # Validate energy matrix
    if E.ndim != 2 or E.shape[0] != E.shape[1]:
        raise ValueError("E must be a square 2D array")
    N_eff = int(N) if N is not None else int(E.shape[0])

    # Normalize solver options
    if solver_opts is None:
        solver_opts_native = _qss.CombinedSolverOptions()
    elif isinstance(solver_opts, _qss.CombinedSolverOptions):
        solver_opts_native = solver_opts
    else:
        solver_opts_native = solver_opts.to_qss()

    records = []
    index_keys = []

    # Loop over all parameter combinations
    for distro in distro_iter:
        native_dist = _coerce_native_distro(distro)
        for beta_val in beta_vals:
            beta_f = float(beta_val)
            for h_val in H_a_vals:
                h_f = float(h_val)
                # Call vectorized backend (loops over n internally)
                results = _qss.chemical_potential_B_multi_n(
                    native_dist,
                    beta_f,
                    n_vals,
                    N_eff,
                    E,
                    h_f,
                    float(g),
                    float(muB_eV_per_T),
                    solver_opts_native,
                    float(a0),
                    float(b0),
                    float(x0),
                )
                for n_val, (mu, iters, conv, err) in zip(n_vals, results):
                    records.append(
                        {
                            "mu": float(mu),
                            "iterations": int(iters),
                            "converged": bool(conv),
                            "residual": float(err),
                        }
                    )
                    index_keys.append((_dist_label(distro), n_val, 1/beta_val, h_val))

    idx = pd.MultiIndex.from_tuples(
        index_keys, names=["distro", "n", "k_BT", "H_a"]
    )
    return pd.DataFrame.from_records(records, index=idx)


chemical_potential_dataframe = chemical_potential

__all__ = ["chemical_potential", "chemical_potential_dataframe"]
