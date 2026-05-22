"""Streamlit widgets for Brent / SA solver configuration."""

from __future__ import annotations

import numpy as np
import streamlit as st

from gui.solver_options import build_brent_options, build_sa_options
from qss.solver_options import CoolingSchedule

SCHED_MAP = {
    "Exponential": CoolingSchedule.Exponential,
    "Linear": CoolingSchedule.Linear,
    "Logarithmic": CoolingSchedule.Logarithmic,
    "Sawtooth": CoolingSchedule.Sawtooth,
}


def solver_controls(container: st.delta_generator.DeltaGenerator) -> tuple[str, object | None, float, float, float]:
    solver = container.selectbox("Algorithm", ["Brent", "SA"], index=0, key="solver_algo")
    opts: object | None = None

    container.markdown("**Search / guess**")
    c1, c2 = container.columns(2)
    search_min = float(c1.number_input("Bracket a0", value=-10.0, key="solver_a0"))
    search_max = float(c2.number_input("Bracket b0", value=10.0, key="solver_b0"))
    use_random = container.checkbox("Random initial guess", value=False, key="solver_rand")
    if use_random:
        x0 = float(np.random.uniform(search_min, search_max))
        container.caption(f"Using x0 = {x0:.6g}")
    else:
        x0 = float(container.number_input("Initial guess x0", value=0.0, key="solver_x0"))

    if solver == "SA":
        container.markdown("**Simulated annealing**")
        max_iters = container.number_input("Max iterations", value=50000, key="sa_max_iters")
        min_iters = container.number_input("Min iterations", value=2000, key="sa_min_iters")
        stall_iters = container.number_input("Stall iterations", value=1000, key="sa_stall")
        T0 = container.number_input("T0", value=1000.0, key="sa_T0")
        T_min = container.number_input("T_min", value=1e-9, format="%.1e", key="sa_Tmin")
        alpha = container.number_input(
            "alpha", value=0.999, min_value=1e-6, max_value=1.0, key="sa_alpha"
        )
        step_decay = container.number_input(
            "step_decay", value=0.999, min_value=1e-6, max_value=1.0, key="sa_step_decay"
        )
        step_sigma = container.number_input("step_sigma", value=5.0, key="sa_sigma")
        sn = container.selectbox("Cooling schedule", list(SCHED_MAP.keys()), key="sa_sched")
        scale_sigma = container.checkbox("Scale sigma by T", value=True, key="sa_scale_sigma")
        reflect_bounds = container.checkbox("Reflect bounds", value=True, key="sa_reflect")
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
        container.markdown("**Brent / bracketing**")
        br1, br2, br3 = container.columns(3)
        abs_tol = float(br1.number_input("abs_tol", value=1e-9, format="%.1e", key="br_abs"))
        rel_tol = float(br2.number_input("rel_tol", value=1e-6, format="%.1e", key="br_rel"))
        max_iters = int(br3.number_input("max_iters", value=5000, key="br_max_iters"))
        br4, br5, br6 = container.columns(3)
        expand_factor = float(br4.number_input("expand_factor", value=2.0, key="br_expand"))
        max_expands = int(br5.number_input("max_expands", value=164, key="br_max_expands"))
        initial_step = float(br6.number_input("initial_step", value=1.0, key="br_initial"))
        brk = container.selectbox("Bracketer", ["Expansion", "Bisection"], key="br_bracketer")
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
