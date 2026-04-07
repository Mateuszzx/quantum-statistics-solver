#!/usr/bin/env python3
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

import qss
from qss.chemical_potential import chemical_potential
from qss.distributions import DistType
from qss.solver_options import brent_defaults


REPO_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class SolverResult:
    mu: float
    residual: float
    iterations: int
    converged: bool


def tight_binding_dispersion(kx: np.ndarray, ky: np.ndarray, t: float, tp: float) -> np.ndarray:
    """Reference 2D tight-binding dispersion used in the comparison scripts."""
    return 2.0 * t * (np.cos(kx) + np.cos(ky)) + 4.0 * tp * np.cos(kx) * np.cos(ky)


def make_reference_energy(N: int, t: float, tp: float) -> np.ndarray:
    """Reference grid used in the older standalone scripts: k in [0, 2pi)."""
    kx = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False, dtype=np.float64)
    ky = np.linspace(0.0, 2.0 * np.pi, N, endpoint=False, dtype=np.float64)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    return tight_binding_dispersion(KX, KY, t, tp)


def make_notebook_energy(N: int, t: float, tp: float) -> np.ndarray:
    """
    Exact grid used in notebooks/Graphs/paramagnetic_mu_appb.ipynb.

    Note: this produces arrays of length N + 1.
    """
    kx = np.linspace(0.0, 2.0 * np.pi * (N - 1) / N, N + 1, dtype=np.float64) - np.pi
    ky = np.linspace(0.0, 2.0 * np.pi * (N - 1) / N, N + 1, dtype=np.float64) - np.pi
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    return tight_binding_dispersion(KX, KY, t, tp)


def emulate_binding_matrix_copy(E: np.ndarray, N: int) -> np.ndarray:
    """
    Reproduce src/bindings/module_chemical_potential.cpp behavior when E.shape != (N, N).

    The binding currently copies the first N*N flat elements and reshapes them as N x N.
    """
    flat = np.ascontiguousarray(E, dtype=np.float64).ravel()
    if flat.size < N * N:
        raise ValueError(f"E does not have enough elements for N={N}: {E.shape}")
    return flat[: N * N].reshape(N, N)


def fd_density(E: np.ndarray, mu: float, T: float) -> float:
    """Total FD density including spin degeneracy."""
    arg = np.clip((E - mu) / T, -700.0, 700.0)
    n_single = 1.0 / (1.0 + np.exp(arg))
    return float(2.0 * np.mean(n_single))


def solve_reference_fd_mu(
    E: np.ndarray,
    T: float,
    n_target: float,
    lo: float = -20.0,
    hi: float = 20.0,
    tol: float = 1e-13,
    max_iter: int = 400,
) -> float:
    """Simple bisection for the monotonic FD density equation."""
    def f(mu: float) -> float:
        return fd_density(E, mu, T) - n_target

    flo = f(lo)
    fhi = f(hi)
    if flo == 0.0:
        return lo
    if fhi == 0.0:
        return hi
    if flo * fhi > 0.0:
        raise RuntimeError(
            f"Reference FD bracket failed: f({lo})={flo:.6e}, f({hi})={fhi:.6e}"
        )

    a = lo
    b = hi
    for _ in range(max_iter):
        mid = 0.5 * (a + b)
        fmid = f(mid)
        if abs(fmid) < tol or abs(b - a) < tol:
            return mid
        if flo * fmid <= 0.0:
            b = mid
            fhi = fmid
        else:
            a = mid
            flo = fmid
    return 0.5 * (a + b)


def solve_qss_fd_mu(E: np.ndarray, T: float, n_target: float, N_arg: int) -> SolverResult:
    """Run the current QSS FD chemical potential solver exactly through the public Python API."""
    py_opts = brent_defaults(
        abs_tol=1e-9,
        rel_tol=1e-6,
        expand_factor=2.0,
        initial_step=1.0,
        max_expands=164,
        max_iters=5000,
    )
    native_opts = py_opts.to_qss()

    df = chemical_potential(
        distros=[DistType.FD],
        n_vals=[n_target],
        kBT_vals=np.array([T], dtype=np.float64),
        H_a_vals=np.array([0.0], dtype=np.float64),
        N=N_arg,
        E=np.ascontiguousarray(E, dtype=np.float64),
        solver_opts=native_opts,
        x0=-1.0,
        a0=-5.0,
        b0=5.0,
    )
    row = df.iloc[0]
    return SolverResult(
        mu=float(row["mu"]),
        residual=float(row["residual"]),
        iterations=int(row["iterations"]),
        converged=bool(row["converged"]),
    )


def print_header() -> None:
    print("FD chemical potential comparison")
    print(f"repo_root: {REPO_ROOT}")
    print(f"qss_python: {qss.__file__}")
    print(f"qss_native: {qss._qss.__file__}")
    print()


def print_case(
    label: str,
    *,
    E: np.ndarray,
    T: float,
    n_target: float,
    qss_result: SolverResult | None = None,
    ref_mu: float | None = None,
    effective_E: np.ndarray | None = None,
) -> None:
    print(label)
    print(f"  E.shape = {E.shape}")
    if ref_mu is not None:
        print(f"  reference_mu = {ref_mu:.12f}")
        print(f"  reference_density(reference_mu) = {fd_density(E, ref_mu, T):.12f}")
    if qss_result is not None:
        print(f"  qss_mu = {qss_result.mu:.12f}")
        print(f"  qss_residual = {qss_result.residual:.12e}")
        print(f"  qss_iterations = {qss_result.iterations}")
        print(f"  qss_converged = {qss_result.converged}")
        print(f"  density_on_visible_E(qss_mu) = {fd_density(E, qss_result.mu, T):.12f}")
        if effective_E is not None:
            print(
                f"  density_on_effective_binding_E(qss_mu) = "
                f"{fd_density(effective_E, qss_result.mu, T):.12f}"
            )
        if ref_mu is not None:
            print(f"  |qss_mu - reference_mu| = {abs(qss_result.mu - ref_mu):.12e}")
    print()


def main() -> None:
    T = 0.01
    n_target = 0.8
    N = 100
    t = -1.0
    tp_reference = 0.25
    tp_notebook = 0.0

    print_header()
    print(
        "Target case: FD, T=0.01, n=0.8, h=0.0, N=100, "
        f"t={t}, tp_reference={tp_reference}"
    )
    print()

    E_reference = make_reference_energy(N, t, tp_reference)
    mu_reference = solve_reference_fd_mu(E_reference, T, n_target)
    qss_reference = solve_qss_fd_mu(E_reference, T, n_target, N_arg=N)
    print_case(
        "Case 1: reference grid, tp=0.25, N=100",
        E=E_reference,
        T=T,
        n_target=n_target,
        qss_result=qss_reference,
        ref_mu=mu_reference,
    )

    E_notebook_tp025 = make_notebook_energy(N, t, tp_reference)
    mu_notebook_tp025 = solve_reference_fd_mu(E_notebook_tp025, T, n_target)
    qss_notebook_tp025 = solve_qss_fd_mu(E_notebook_tp025, T, n_target, N_arg=E_notebook_tp025.shape[0])
    print_case(
        "Case 2: notebook grid, tp=0.25, N passed correctly as 101",
        E=E_notebook_tp025,
        T=T,
        n_target=n_target,
        qss_result=qss_notebook_tp025,
        ref_mu=mu_notebook_tp025,
    )

    E_notebook_tp0 = make_notebook_energy(N, t, tp_notebook)
    mu_notebook_tp0 = solve_reference_fd_mu(E_notebook_tp0, T, n_target)
    qss_notebook_tp0 = solve_qss_fd_mu(E_notebook_tp0, T, n_target, N_arg=E_notebook_tp0.shape[0])
    print_case(
        "Case 3: notebook grid, tp=0.0, N passed correctly as 101",
        E=E_notebook_tp0,
        T=T,
        n_target=n_target,
        qss_result=qss_notebook_tp0,
        ref_mu=mu_notebook_tp0,
    )

    E_effective_bug = emulate_binding_matrix_copy(E_notebook_tp0, N)
    mu_effective_bug = solve_reference_fd_mu(E_effective_bug, T, n_target)
    qss_notebook_exact = solve_qss_fd_mu(E_notebook_tp0, T, n_target, N_arg=N)
    print_case(
        "Case 4: exact notebook call, tp=0.0, E is 101x101 but N=100",
        E=E_notebook_tp0,
        T=T,
        n_target=n_target,
        qss_result=qss_notebook_exact,
        ref_mu=mu_effective_bug,
        effective_E=E_effective_bug,
    )


if __name__ == "__main__":
    main()
