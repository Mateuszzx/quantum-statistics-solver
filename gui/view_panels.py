"""Panels view: parameter sweeps, experiment run, sweep plots."""

from __future__ import annotations

import traceback

import numpy as np
import streamlit as st

from gui.constants import SWEEP_VAR_LABELS
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
    n_vals, kBT_vals, H_a_vals, U_vals = experiment_parameter_sweeps(st.sidebar)

    with st.sidebar.expander("Solver parameters", expanded=False) as solver_exp:
        solver, solver_opts, a0, b0, x0 = solver_controls(solver_exp)

    logs_dir = SOLVER_ROOT / "logs"
    out_csv = SOLVER_ROOT / "data" / "gui_results.csv"
    logs_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    est = len(distros) * len(n_vals) * len(kBT_vals) * len(H_a_vals) * len(U_vals)
    st.sidebar.metric("Total solves (approx.)", f"{est:,}")

    if st.sidebar.button("Run experiment", type="primary"):
        if not distros:
            st.error("Select at least one distribution.")
            return
        prog_cb, prog_ph = st_progress_pair(st)
        try:
            e_k, _kx, _ky = build_dispersion(dim, n_grid, t, t_prim)
            with st.spinner("Running chemical potential solves…"):
                mu_df = run_experiment(
                    n_vals=n_vals,
                    kBT_vals=kBT_vals,
                    h_a_vals=H_a_vals,
                    distros=distros,
                    e_k=e_k,
                    N=int(e_k.shape[0]),
                    solver=solver,
                    solver_opts=solver_opts,
                    logs_dir=str(logs_dir),
                    output_file=str(out_csv),
                    x0=x0,
                    a0=a0,
                    b0=b0,
                    U_vals=U_vals,
                    progress=prog_cb,
                )
            prog_ph.empty()
            if mu_df.empty:
                st.error("No results returned. Check stderr / logs for solver errors.")
                return

            st.session_state["meta"] = {
                "n_vals": np.array(n_vals, dtype=float),
                "kBT_vals": np.array(kBT_vals, dtype=float),
                "H_a_vals": np.array(H_a_vals, dtype=float),
                "U_vals": np.array(U_vals, dtype=float),
                "distros": distros,
                "t_prim": t_prim,
                "N": int(e_k.shape[0]),
                "dim": dim,
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

    st.header("Sweep panels")
    q_name = st.selectbox(
        "Quantity", options=derived_display_names(), index=0, key="panel_qty"
    )
    x_var = st.selectbox(
        "X axis variable", options=["n", "k_BT", "H_a", "U"], index=0, key="panel_x"
    )
    x_scale = st.radio(
        "X scale", ["linear", "log"], horizontal=True, index=0, key="panel_xscale"
    )
    panel_var_raw = st.selectbox(
        "Sweep in panels",
        options=["(none)", "n", "k_BT", "H_a", "U"],
        index=2,
        key="panel_pan",
    )
    panel_var = None if panel_var_raw == "(none)" else str(panel_var_raw)
    if panel_var == x_var:
        st.warning("Panel variable must be different from X axis. Switching panels to (none).")
        panel_var = None

    n_cols = int(style["cols"])

    fixed: dict[str, float] = {}
    for var, opts in (("n", n_arr), ("k_BT", kBT_arr), ("H_a", H_arr), ("U", U_arr)):
        if var in (x_var, panel_var):
            continue
        vals = list(np.array(opts, dtype=float))
        if not vals:
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
