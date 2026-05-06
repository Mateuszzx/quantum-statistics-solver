"""Map Streamlit / UI choices to ``CombinedSolverOptions`` for ``qss``."""

from __future__ import annotations

from qss.solver_options import BracketerType, CombinedSolverOptions, CoolingSchedule, brent_defaults, sa_defaults


def build_sa_options(
    *,
    max_iters: int,
    min_iters: int,
    stall_iters: int,
    T0: float,
    T_min: float,
    alpha: float,
    step_decay: float,
    step_sigma: float,
    sched: CoolingSchedule,
    search_min: float,
    search_max: float,
    scale_sigma_by_t: bool,
    reflect_bounds: bool,
) -> CombinedSolverOptions:
    """Return the Python-side options object (``chemical_potential`` calls ``to_qss()`` internally)."""
    py = sa_defaults(
        T0=float(T0),
        T_min=float(T_min),
        alpha=float(alpha),
        max_iters=int(max_iters),
        min_iters=int(min_iters),
        stall_iters=int(stall_iters),
        step_decay=float(step_decay),
        step_sigma=float(step_sigma),
        bounds=[(float(search_min), float(search_max))],
        scale_sigma_by_T=bool(scale_sigma_by_t),
        reflect_bounds=bool(reflect_bounds),
        sched=sched,
    )
    return py


def build_brent_options(
    *,
    abs_tol: float,
    rel_tol: float,
    max_iters: int,
    expand_factor: float,
    max_expands: int,
    initial_step: float,
    bracketer: str,
) -> CombinedSolverOptions:
    """Return the Python-side options object."""
    br = BracketerType.Expansion if bracketer == "Expansion" else BracketerType.Bisection
    py = brent_defaults(
        abs_tol=float(abs_tol),
        rel_tol=float(rel_tol),
        max_iters=int(max_iters),
        bracketer=br,
        expand_factor=float(expand_factor),
        max_expands=int(max_expands),
        initial_step=float(initial_step),
    )
    return py

