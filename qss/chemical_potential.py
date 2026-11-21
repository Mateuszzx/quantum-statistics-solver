"""High-level wrappers for chemical potential calculations."""

from __future__ import annotations

from typing import Sequence

import numpy as np
import pandas as pd

try:
    from . import _qss  # type: ignore[attr-defined, unused-ignore]
except ImportError as e:
    raise ImportError("The _qss extension module is not available.") from e

from .distributions import DistType as PyDistType
from .distributions import distribution


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


def band_filling_from_mu(  # noqa: PLR0913
    E: np.ndarray,
    *,
    mu: float,
    kBT: float,
    H_a: float,
    distro: PyDistType | _qss.DistType | str,
    U: float = 0.0,
    g: float = 2.0,
    muB_eV_per_T: float = 1.0,
) -> float:
    """Band filling ``n`` implied by a prescribed chemical potential (same as C++ residual)."""
    delta = 0.5 * float(g) * float(muB_eV_per_T) * float(H_a)
    label = _dist_label(distro)
    n_up = distribution(E, float(mu), float(kBT), delta, 1, label, U=float(U))
    n_dn = distribution(E, float(mu), float(kBT), delta, -1, label, U=float(U))
    n_up = np.clip(np.asarray(n_up, dtype=float), 0.0, 1.0)
    n_dn = np.clip(np.asarray(n_dn, dtype=float), 0.0, 1.0)
    total = n_up + n_dn
    if label not in ("FD", "FLB", "HUBBARD"):
        total = np.clip(total, 0.0, 1.0)
    return float(np.mean(np.abs(total)))


def chemical_potential(  # noqa: PLR0913, PLR0912
    distros: Sequence[PyDistType | _qss.DistType | str],
    n_vals: Sequence[float] | None = None,
    kBT_vals: Sequence[float] = (),
    H_a_vals: Sequence[float] = (),
    E: np.ndarray | None = None,
    U_vals: Sequence[float] | float = 0.0,
    N: int | None = None,
    solver_opts: _qss.CombinedSolverOptions | None = None,
    g: float = 2.0,
    muB_eV_per_T: float = 1.0,
    a0: float = -10.0,
    b0: float = 10.0,
    x0: float = 0.0,
    mu_vals: Sequence[float] | None = None,
) -> pd.DataFrame:
    """
    Compute chemical potential for multiple parameter combinations.

    Parameters
    ----------
    distros : distribution type or sequence of distribution types
        FD, FLB, NFLB, or ND (as string, enum, or native type)
    n_vals : sequence of float, optional
        Band filling values (solve μ from *n*). Required unless *mu_vals* is set.
    kBT_vals : sequence of float
        Temperature values in energy units (k_B * T)
    H_a_vals : sequence of float
        Applied magnetic field values (Tesla)
    E : ndarray
        Energy matrix (square 2D array)
    U_vals : sequence of float or float
        Hubbard U parameter(s) (default: 0.0)
    N : int, optional
        Effective system size (defaults to E.shape[0])
    solver_opts : CombinedSolverOptions, optional
        Solver configuration (ignored when *mu_vals* is set)
    g : float
        Landé g-factor (default: 2.0)
    muB_eV_per_T : float
        Bohr magneton in eV/T (default: 5.788381806e-5)
    a0, b0, x0 : float
        Initial bracketing and guess parameters
    mu_vals : sequence of float, optional
        Prescribe μ explicitly and skip the μ(n) root solve. Band filling *n* is
        computed from the occupancy for each point (same convention as the solver residual).

    Returns
    -------
    pd.DataFrame
        MultiIndex DataFrame with columns: mu, iterations, converged, residual
        Index levels: distro, n, kBT, H_a
    """
    if E is None:
        raise ValueError("E (energy matrix) is required")
    if isinstance(distros, (str, _qss.DistType, PyDistType)):
        distro_iter = [distros]
    else:
        distro_iter = list(distros)

    U_iter = [float(U_vals)] if isinstance(U_vals, (int, float)) else list(U_vals)

    # Validate energy matrix
    if E.ndim != 2 or E.shape[0] != E.shape[1]:
        raise ValueError("E must be a square 2D array")
    N_eff = int(N) if N is not None else int(E.shape[0])
    if E.shape != (N_eff, N_eff):
        raise ValueError(f"E shape {E.shape} does not match N={N_eff}; expected E.shape == (N, N).")

    prescribe_mu = mu_vals is not None
    if prescribe_mu and n_vals is not None:
        raise ValueError("Pass either n_vals (solve μ) or mu_vals (prescribe μ), not both.")
    if not prescribe_mu and n_vals is None:
        raise ValueError("Either n_vals or mu_vals must be provided.")

    records = []
    index_keys = []

    if prescribe_mu:
        for distro in distro_iter:
            for kBT_val in kBT_vals:
                kBT_f = float(kBT_val)
                for h_val in H_a_vals:
                    h_f = float(h_val)
                    for u_val in U_iter:
                        u_f = float(u_val)
                        for mu_val in mu_vals:  # type: ignore[union-attr]
                            mu_f = float(mu_val)
                            n_f = band_filling_from_mu(
                                E,
                                mu=mu_f,
                                kBT=kBT_f,
                                H_a=h_f,
                                distro=distro,
                                U=u_f,
                                g=float(g),
                                muB_eV_per_T=float(muB_eV_per_T),
                            )
                            records.append(
                                {
                                    "mu": mu_f,
                                    "iterations": 0,
                                    "converged": True,
                                    "residual": 0.0,
                                }
                            )
                            index_keys.append((_dist_label(distro), n_f, kBT_val, h_val, u_val))
        idx = pd.MultiIndex.from_tuples(index_keys, names=["distro", "n", "k_BT", "H_a", "U"])
        return pd.DataFrame.from_records(records, index=idx)

    # Normalize solver options
    if solver_opts is None:
        solver_opts_native = _qss.CombinedSolverOptions()
    elif isinstance(solver_opts, _qss.CombinedSolverOptions):
        solver_opts_native = solver_opts
    elif isinstance(solver_opts, _qss.SimulatedAnnealingOptions):
        solver_opts_native = _qss.CombinedSolverOptions()
        solver_opts_native.solver = _qss.SolverType.SimulatedAnnealing
        solver_opts_native.sa = solver_opts
    else:
        # Try to convert from Python object
        native = solver_opts.to_qss()
        if isinstance(native, _qss.SimulatedAnnealingOptions):
            solver_opts_native = _qss.CombinedSolverOptions()
            solver_opts_native.solver = _qss.SolverType.SimulatedAnnealing
            solver_opts_native.sa = native
        elif isinstance(native, _qss.CombinedSolverOptions):
            solver_opts_native = native
        else:
            # Assume RootSolverOptions
            solver_opts_native = _qss.CombinedSolverOptions()
            solver_opts_native.root = native

    # Loop over all parameter combinations
    for distro in distro_iter:
        native_dist = _coerce_native_distro(distro)
        for kBT_val in kBT_vals:
            kBT_f = float(kBT_val)
            for h_val in H_a_vals:
                h_f = float(h_val)
                for u_val in U_iter:
                    u_f = float(u_val)
                    # Call vectorized backend (loops over n internally)
                    results = _qss.chemical_potential_B_multi_n(
                        native_dist,
                        kBT_f,
                        n_vals,
                        N_eff,
                        E,
                        h_f,
                        u_f,
                        float(g),
                        float(muB_eV_per_T),
                        solver_opts_native,
                        float(a0),
                        float(b0),
                        float(x0),
                    )
                    for n_val, (mu, iters, conv, err) in zip(n_vals, results):  # type: ignore[arg-type]
                        records.append(
                            {
                                "mu": float(mu),
                                "iterations": int(iters),
                                "converged": bool(conv),
                                "residual": float(err),
                            }
                        )
                        index_keys.append((_dist_label(distro), n_val, kBT_val, h_val, u_val))

    idx = pd.MultiIndex.from_tuples(index_keys, names=["distro", "n", "k_BT", "H_a", "U"])
    return pd.DataFrame.from_records(records, index=idx)


chemical_potential_dataframe = chemical_potential

__all__ = ["band_filling_from_mu", "chemical_potential", "chemical_potential_dataframe"]
