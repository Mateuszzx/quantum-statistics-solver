"""Streamlit widgets for Brent / SA solver configuration."""

from __future__ import annotations

import numpy as np
import streamlit as st

from gui.defaults import cfg, option_index
from gui.solver_options import build_brent_options, build_sa_options
from qss.solver_options import CoolingSchedule

SCHED_MAP = {
    "Exponential": CoolingSchedule.Exponential,
    "Linear": CoolingSchedule.Linear,
    "Logarithmic": CoolingSchedule.Logarithmic,
    "Sawtooth": CoolingSchedule.Sawtooth,
}

_BRACKETERS = ["Expansion", "Bisection"]


def solver_controls(container: st.delta_generator.DeltaGenerator) -> tuple[str, object | None, float, float, float]:
    solver = container.selectbox(
        "Algorithm",
        ["Brent", "SA"],
        index=int(cfg("solver", "algorithm_index")),
        key="solver_algo",
    )
    opts: object | None = None
    search = cfg("solver", "search")

    container.markdown("**Search / guess**")
    c1, c2 = container.columns(2)
    search_min = float(c1.number_input("Bracket a0", value=float(search["a0"]), key="solver_a0"))
    search_max = float(c2.number_input("Bracket b0", value=float(search["b0"]), key="solver_b0"))
    use_random = container.checkbox(
        "Random initial guess", value=bool(search["random_guess"]), key="solver_rand"
    )
    if use_random:
        x0 = float(np.random.uniform(search_min, search_max))
        container.caption(f"Using x0 = {x0:.6g}")
    else:
        x0 = float(container.number_input("Initial guess x0", value=float(search["x0"]), key="solver_x0"))

    if solver == "SA":
        sa = cfg("solver", "sa")
        container.markdown("**Simulated annealing**")
        max_iters = container.number_input("Max iterations", value=int(sa["max_iters"]), key="sa_max_iters")
        min_iters = container.number_input("Min iterations", value=int(sa["min_iters"]), key="sa_min_iters")
        stall_iters = container.number_input("Stall iterations", value=int(sa["stall_iters"]), key="sa_stall")
        T0 = container.number_input("T0", value=float(sa["T0"]), key="sa_T0")
        T_min = container.number_input("T_min", value=float(sa["T_min"]), format="%.1e", key="sa_Tmin")
        alpha = container.number_input(
            "alpha", value=float(sa["alpha"]), min_value=1e-6, max_value=1.0, key="sa_alpha"
        )
        step_decay = container.number_input(
            "step_decay", value=float(sa["step_decay"]), min_value=1e-6, max_value=1.0, key="sa_step_decay"
        )
        step_sigma = container.number_input("step_sigma", value=float(sa["step_sigma"]), key="sa_sigma")
        sched_names = list(SCHED_MAP.keys())
        sn = container.selectbox(
            "Cooling schedule",
            sched_names,
            index=option_index(sched_names, str(sa["schedule"]), field="solver.sa.schedule"),
            key="sa_sched",
        )
        scale_sigma = container.checkbox(
            "Scale sigma by T", value=bool(sa["scale_sigma_by_t"]), key="sa_scale_sigma"
        )
        reflect_bounds = container.checkbox(
            "Reflect bounds", value=bool(sa["reflect_bounds"]), key="sa_reflect"
        )
        opts = build_sa_options(
            max_iters=int(max_iters),
            min_iters=int(min_iters),
            stall_iters=int(stall_iters),
            T0=float(T0),
            T_min=float(T_min),
            alpha=float(alpha),
            step_decay=float(step_decay),
            step_sigma=float(step_sigma),
            sched=SCHED_MAP[str(sn)],
            search_min=float(search_min),
            search_max=float(search_max),
            scale_sigma_by_t=bool(scale_sigma),
            reflect_bounds=bool(reflect_bounds),
        )
    else:
        br = cfg("solver", "brent")
        container.markdown("**Brent / bracketing**")
        br1, br2, br3 = container.columns(3)
        abs_tol = float(br1.number_input("abs_tol", value=float(br["abs_tol"]), format="%.1e", key="br_abs"))
        rel_tol = float(br2.number_input("rel_tol", value=float(br["rel_tol"]), format="%.1e", key="br_rel"))
        max_iters = int(br3.number_input("max_iters", value=int(br["max_iters"]), key="br_max_iters"))
        br4, br5, br6 = container.columns(3)
        expand_factor = float(br4.number_input("expand_factor", value=float(br["expand_factor"]), key="br_expand"))
        max_expands = int(br5.number_input("max_expands", value=int(br["max_expands"]), key="br_max_expands"))
        initial_step = float(br6.number_input("initial_step", value=float(br["initial_step"]), key="br_initial"))
        brk = container.selectbox(
            "Bracketer",
            _BRACKETERS,
            index=option_index(_BRACKETERS, str(br["bracketer"]), field="solver.brent.bracketer"),
            key="br_bracketer",
        )
        opts = build_brent_options(
            abs_tol=float(abs_tol),
            rel_tol=float(rel_tol),
            max_iters=int(max_iters),
            expand_factor=float(expand_factor),
            max_expands=int(max_expands),
            initial_step=float(initial_step),
            bracketer=str(brk),
        )

    return solver, opts, search_min, search_max, x0
