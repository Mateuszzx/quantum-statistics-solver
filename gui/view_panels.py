"""Panels view: parameter sweeps, experiment run, sweep plots."""

from __future__ import annotations

import traceback

import numpy as np
import streamlit as st

from gui.constants import SWEEP_VAR_LABELS
from gui.defaults import cfg
from gui.derived_quantities import (
    derived_compute_summary,
    derived_display_names,
    derived_quantities_ready,
    df_for_quantity,
    store_derived_in_session,
)
from gui.lattice import build_dispersion
from gui.panel_plot import panel_plot
from gui.paths import quantum_statistics_solver_root
from gui.plot_controls import plot_style_controls
from gui.sidebar_controls import sidebar_distros, sidebar_lattice
from gui.solver_controls import solver_controls
from gui.streamlit_utils import show_current_figure, st_progress_pair
from gui.sweep_input import experiment_parameter_sweeps
from utils.experiment_utils import run_experiment

SOLVER_ROOT = quantum_statistics_solver_root()


def view_panels() -> None:  # noqa: PLR0912, PLR0915
    distros = sidebar_distros(key_prefix="panel_")
    dim, n_grid, t, t_prim = sidebar_lattice(key_prefix="panel_")

    st.sidebar.subheader("Experiment parameters")
    prescribe_mu = st.sidebar.checkbox(
        r"Prescribe $\mu$ (skip $\mu(n)$ solve)",
        value=bool(cfg("sweeps", "prescribe_mu")),
        key="panel_prescribe_mu",
        help="Use an explicit chemical potential sweep instead of solving μ from band filling n.",
    )
    sweeps = experiment_parameter_sweeps(st.sidebar, prescribe_mu=prescribe_mu)
    kBT_vals = sweeps["k_BT"]
    H_a_vals = sweeps["H_a"]
    U_vals = sweeps["U"]
    n_vals = sweeps.get("n")
    mu_vals = sweeps.get("mu")

    solver = "Brent"
    solver_opts = a0 = b0 = x0 = None
    if prescribe_mu:
        st.sidebar.caption("Solver parameters are unused when μ is prescribed.")
    else:
        with st.sidebar.expander("Solver parameters", expanded=False) as solver_exp:
            solver, solver_opts, a0, b0, x0 = solver_controls(solver_exp)

    logs_dir = SOLVER_ROOT / "logs"
    out_csv = SOLVER_ROOT / "data" / "gui_results.csv"
    logs_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    primary_len = len(mu_vals) if prescribe_mu else len(n_vals)  # type: ignore[arg-type]
    est = len(distros) * primary_len * len(kBT_vals) * len(H_a_vals) * len(U_vals)
    st.sidebar.metric("Total solves (approx.)", f"{est:,}")

    if st.sidebar.button("Run experiment", type="primary"):
        if not distros:
            st.error("Select at least one distribution.")
            return
        prog_cb, prog_ph = st_progress_pair(st)
        try:
            e_k, _kx, _ky = build_dispersion(dim, n_grid, t, t_prim)
            spinner_msg = (
                "Evaluating occupancy at prescribed μ…"
                if prescribe_mu
                else "Running chemical potential solves…"
            )
            with st.spinner(spinner_msg):
                run_kwargs: dict[str, object] = {
                    "kBT_vals": kBT_vals,
                    "h_a_vals": H_a_vals,
                    "distros": distros,
                    "e_k": e_k,
                    "N": int(e_k.shape[0]),
                    "logs_dir": str(logs_dir),
                    "output_file": str(out_csv),
                    "U_vals": U_vals,
                    "progress": prog_cb,
                }
                if prescribe_mu:
                    run_kwargs["mu_vals"] = mu_vals
                else:
                    run_kwargs.update(
                        {
                            "n_vals": n_vals,
                            "solver": solver,
                            "solver_opts": solver_opts,
                            "x0": x0,
                            "a0": a0,
                            "b0": b0,
                        }
                    )
                mu_df = run_experiment(**run_kwargs)  # type: ignore[arg-type]
            prog_ph.empty()
            if mu_df.empty:
                st.error("No results returned. Check stderr / logs for solver errors.")
                return

            n_from_df = (
                np.array(sorted(mu_df.index.get_level_values("n").unique()), dtype=float)
                if "n" in mu_df.index.names
                else np.array([], dtype=float)
            )
            st.session_state["meta"] = {
                "n_vals": n_from_df if prescribe_mu else np.array(n_vals, dtype=float),
                "mu_vals": np.array(mu_vals, dtype=float) if prescribe_mu else None,
                "kBT_vals": np.array(kBT_vals, dtype=float),
                "H_a_vals": np.array(H_a_vals, dtype=float),
                "U_vals": np.array(U_vals, dtype=float),
                "distros": distros,
                "t_prim": t_prim,
                "N": int(e_k.shape[0]),
                "dim": dim,
                "prescribe_mu": prescribe_mu,
            }
            with st.spinner(f"Computing {derived_compute_summary()}…"):
                store_derived_in_session(st.session_state, mu_df=mu_df, e_k=e_k)
            st.success("Finished.")
        except Exception as e:
            prog_ph.empty()
            st.error(str(e))
            st.code(traceback.format_exc())
            return

    if st.session_state.get("mu_df") is None:
        st.info("Configure the sidebar and click **Run experiment** to plot sweep panels.")
        return

    if not derived_quantities_ready(st.session_state):
        st.error("Derived quantities missing — run the experiment again.")
        return

    with st.expander("Plot style", expanded=False):
        style = plot_style_controls()

    meta = st.session_state["meta"]
    distros = meta["distros"]
    n_arr = meta["n_vals"]
    kBT_arr = meta["kBT_vals"]
    H_arr = meta["H_a_vals"]
    U_arr = meta["U_vals"]
    mu_arr = meta.get("mu_vals")
    used_prescribe = bool(meta.get("prescribe_mu"))

    st.header("Sweep panels")
    mu_count = int(len(np.asarray(mu_arr))) if mu_arr is not None else 0
    if used_prescribe:
        st.caption(
            r"Prescribed $\mu$: occupancy-implied $n$ is available as an axis, "
            r"but is **not** used as a Fixed filter (it changes with $T$, $H_a$, $U$). "
            r"For **constant-$\mu$ panels**, set $\mu$ to Scalar in the sidebar, "
            r"then e.g. X = $k_B T$, panels = $H_a$."
        )
    p_cfg = cfg("panels")
    # Prefer plotting vs T / H / U when μ is held fixed (Scalar); include μ when swept.
    if used_prescribe:
        sweep_vars = (
            ["mu", "n", "k_BT", "H_a", "U"] if mu_count > 1 else ["k_BT", "H_a", "U", "n", "mu"]
        )
        default_x = 0  # k_BT when single μ, else mu
        default_pan = min(2, len(sweep_vars))  # H_a-ish when single μ; else k_BT
        # Map default_pan into panel_opts which starts with "(none)"
        default_pan_opt = default_pan + 1
    else:
        sweep_vars = ["n", "k_BT", "H_a", "U"]
        default_x = int(p_cfg["x_var_index"])
        default_pan_opt = int(p_cfg["panel_var_index"])
    panel_opts = ["(none)", *sweep_vars]
    q_name = st.selectbox(
        "Quantity",
        options=derived_display_names(),
        index=int(p_cfg["quantity_index"]),
        key="panel_qty",
    )
    x_var = st.selectbox(
        "X axis variable",
        options=sweep_vars,
        index=min(default_x, len(sweep_vars) - 1),
        key="panel_x",
    )
    x_scale = st.radio(
        "X scale",
        ["linear", "log"],
        horizontal=True,
        index=int(p_cfg["x_scale_index"]),
        key="panel_xscale",
    )
    panel_var_raw = st.selectbox(
        "Sweep in panels",
        options=panel_opts,
        index=min(default_pan_opt, len(panel_opts) - 1),
        key="panel_pan",
    )
    panel_var = None if panel_var_raw == "(none)" else str(panel_var_raw)
    if panel_var == x_var:
        st.warning("Panel variable must be different from X axis. Switching panels to (none).")
        panel_var = None

    n_cols = int(style["cols"])

    # Independent controls only — when μ is prescribed, n is a dependent output.
    fixed_sources: list[tuple[str, object]] = [
        ("k_BT", kBT_arr),
        ("H_a", H_arr),
        ("U", U_arr),
    ]
    if used_prescribe and mu_arr is not None:
        fixed_sources.insert(0, ("mu", mu_arr))
    else:
        fixed_sources.insert(0, ("n", n_arr))

    fixed: dict[str, float] = {}
    for var, opts in fixed_sources:
        if var in (x_var, panel_var):
            continue
        vals = list(np.array(opts, dtype=float))
        if not vals:
            continue
        # Single-valued controls need no widget (always that value).
        if len(vals) == 1:
            fixed[var] = float(vals[0])
            continue
        default_val = float(vals[0])
        fixed[var] = float(
            st.selectbox(
                f"Fixed ${var}$",
                options=vals,
                index=int(np.argmin(np.abs(np.array(vals, dtype=float) - default_val))),
                format_func=lambda x: f"{float(x):.6g}",
                key=f"panel_fix_{var}",
            )
        )

    dfq, y_col, y_label = df_for_quantity(st.session_state, display_name=str(q_name))
    panel_plot(
        dfq,
        y_col=y_col,
        y_label=y_label,
        x_var=str(x_var),
        x_scale=str(x_scale),
        panel_var=panel_var,
        n_cols=n_cols,
        fixed=fixed,
        distros=distros,
        style=style,
    )
    show_current_figure(st)
    fixed_md = (
        ", ".join([f"{SWEEP_VAR_LABELS.get(k, k)}={v:.6g}" for k, v in fixed.items()])
        or "(none)"
    )
    x_md = SWEEP_VAR_LABELS.get(str(x_var), str(x_var))
    pan_md = SWEEP_VAR_LABELS.get(str(panel_var_raw), "none") if panel_var_raw != "(none)" else "none"
    st.markdown(
        f"**Quantity**: {q_name}  \n"
        f"**X**: {x_md} ({x_scale})  \n"
        f"**Panels**: {pan_md}  \n"
        f"**Fixed**: {fixed_md}"
    )

    with st.expander("Data (table + download)", expanded=False):
        st.dataframe(dfq.reset_index())
        st.download_button(
            f"Download {y_col} CSV",
            dfq.reset_index().to_csv(index=False).encode("utf-8"),
            f"gui_{y_col}.csv",
            "text/csv",
        )
