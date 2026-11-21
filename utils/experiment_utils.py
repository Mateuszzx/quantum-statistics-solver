"""Batch chemical-potential experiments (used by notebooks and the Streamlit GUI)."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Protocol
import uuid

from numpy.typing import NDArray
import pandas as pd
from tqdm.auto import tqdm

from qss.chemical_potential import chemical_potential
from qss.distributions import DistType
from qss.solver_options import CombinedSolverOptions, CoolingSchedule, brent_defaults, sa_defaults


class _HasLogFile(Protocol):
    log_file: str


class _HasSA(Protocol):
    sa: _HasLogFile


SolverOptsInput = CombinedSolverOptions | _HasLogFile | _HasSA


def _set_log_path(native_opts: SolverOptsInput, path: str) -> None:
    if hasattr(native_opts, "log_file"):
        native_opts.log_file = path
    elif hasattr(native_opts, "sa"):
        native_opts.sa.log_file = path


def _resolve_solver_opts(solver: str, solver_opts: SolverOptsInput | None) -> SolverOptsInput:
    if solver_opts is not None:
        return solver_opts
    s = solver.upper()
    if s == "SA":
        py = sa_defaults(
            T0=1e3,
            T_min=1e-9,
            alpha=0.999,
            max_iters=50000,
            min_iters=2000,
            stall_iters=1000,
            step_decay=0.999,
            step_sigma=5.0,
            bounds=[(-10.0, 10.0)],
            scale_sigma_by_T=True,
            reflect_bounds=True,
            sched=CoolingSchedule.Exponential,
        )
        return py.to_qss()
    if s in ("BRENT", "BRENTQ", "ROOT"):
        return brent_defaults().to_qss()
    raise ValueError(f"Unknown solver: {solver}")


def run_experiment(  # noqa: PLR0912, PLR0913
    n_vals: Sequence[float] | None = None,
    kBT_vals: Sequence[float] = (),
    h_a_vals: Sequence[float] = (),
    distros: Sequence[DistType] = (),
    e_k: NDArray[object] | None = None,
    N: int = 100,
    solver: str = "SA",
    solver_opts: SolverOptsInput | None = None,
    logs_dir: str | Path = "logs",
    output_file: str | Path = "chemical_potential_experiment_results.csv",
    x0: float = 0.0,
    a0: float = -10.0,
    b0: float = 10.0,
    U_vals: Sequence[float] | float = 0.0,
    progress: Callable[[int, int], None] | None = None,
    mu_vals: Sequence[float] | None = None,
) -> pd.DataFrame:
    """
    Run ``chemical_potential`` for a Cartesian grid of (distro, k_B T, H_a, u, n|μ).

    Pass either *n_vals* (solve μ from band filling) or *mu_vals* (prescribe μ and
    skip the root solve). When *mu_vals* is set, solver options are unused.
    """
    if e_k is None:
        raise ValueError("e_k is required")
    if mu_vals is not None and n_vals is not None:
        raise ValueError("Pass either n_vals or mu_vals, not both.")
    if mu_vals is None and n_vals is None:
        raise ValueError("Either n_vals or mu_vals must be provided.")

    prescribe_mu = mu_vals is not None
    native_opts = None if prescribe_mu else _resolve_solver_opts(solver, solver_opts)

    logs_path = Path(logs_dir)
    logs_path.mkdir(parents=True, exist_ok=True)

    U_iter = [float(U_vals)] if isinstance(U_vals, int | float) else list(U_vals)
    sweep_vals = list(mu_vals) if prescribe_mu else list(n_vals)  # type: ignore[arg-type]
    results: list[pd.DataFrame] = []
    total = len(distros) * len(kBT_vals) * len(h_a_vals) * len(sweep_vals) * len(U_iter)
    done = 0

    with tqdm(total=total, desc="Simulations", disable=progress is not None) as pbar:
        for d in distros:
            for kbt in kBT_vals:
                for h in h_a_vals:
                    for u_val in U_iter:
                        for sweep_val in sweep_vals:
                            try:
                                if prescribe_mu:
                                    res_df = chemical_potential(
                                        distros=[d],
                                        mu_vals=[float(sweep_val)],
                                        kBT_vals=[float(kbt)],
                                        H_a_vals=[float(h)],
                                        N=int(N),
                                        E=e_k,
                                        U_vals=[float(u_val)],
                                    )
                                else:
                                    rnd = uuid.uuid4().hex[:6]
                                    dist_name = str(d).split(".")[-1]
                                    log_fn = (
                                        logs_path
                                        / f"log_{solver.lower()}_{dist_name}_kbt{kbt:.4g}_n{sweep_val:.4g}_H{h:.4g}_U{u_val:.4g}_{rnd}.csv"
                                    )
                                    assert native_opts is not None
                                    _set_log_path(native_opts, str(log_fn))
                                    res_df = chemical_potential(
                                        distros=[d],
                                        n_vals=[float(sweep_val)],
                                        kBT_vals=[float(kbt)],
                                        H_a_vals=[float(h)],
                                        N=int(N),
                                        E=e_k,
                                        U_vals=[float(u_val)],
                                        solver_opts=native_opts,
                                        x0=float(x0),
                                        a0=float(a0),
                                        b0=float(b0),
                                    )
                                results.append(res_df)
                            except Exception as e:
                                tag = "mu" if prescribe_mu else "n"
                                print(
                                    f"\nError for: {d}, {tag}={sweep_val}, kBT={kbt}, H={h}, U={u_val}\nDetails: {e}"
                                )
                            done += 1
                            pbar.update(1)
                            if progress is not None:
                                progress(done, total)

    if not results:
        idx = pd.MultiIndex.from_tuples([], names=["distro", "n", "k_BT", "H_a", "U"])
        return pd.DataFrame(index=idx)

    final_df = pd.concat(results)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_file)
    return final_df
