"""
Streamlit GUI for quantum-statistics-solver: reproduces notebook workflows (sweeps in
``n``, ``k_B T``, ``H_a``; μ, entropy, energy, free energy, magnetization, ``n_{kd}``;
distribution / BZ cuts when ``qss-results/tools`` is available).

Run from the ``quantum-statistics-solver`` directory::

    streamlit run gui/web_app.py
"""

from __future__ import annotations

import traceback

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from gui.lattice import build_dispersion
from gui.panel_plot import panel_plot
from gui.paths import quantum_statistics_solver_root
from gui.scientific_plot import ScientificPlotter
from gui.solver_options import build_brent_options, build_sa_options
from gui.streamlit_utils import show_current_figure, st_progress_pair
from gui.sweep_grids import parse_float_list
from gui.sweep_input import experiment_parameter_sweeps
from qss.chemical_potential import chemical_potential
from qss.distributions import DistType, distribution
from qss.energy import energy as compute_energy
from qss.entropy import entropy as compute_entropy
from qss.heat_capacity import heat_capacity as compute_heat_capacity
from qss.magnetization import calculate_magnetization
from qss.nkd import n_kd as compute_nkd
from qss import solver_options as so
from qss.solver_options import CoolingSchedule
from utils.experiment_utils import run_experiment

st.set_page_config(page_title="Quantum Statistics Solver", layout="wide")

SOLVER_ROOT = quantum_statistics_solver_root()

DIM = "2D"
DIM_2D = 2

DIST_MAP: dict[str, DistType] = {
    "FD": DistType.FD,
    "FLB": DistType.FLB,
    "NFL": DistType.NFL,
    "NFLB": DistType.NFLB,
    "SSL": DistType.SSL,
    "HUBBARD": DistType.HUBBARD,
}

SCHED_MAP = {
    "Exponential": CoolingSchedule.Exponential,
    "Linear": CoolingSchedule.Linear,
    "Logarithmic": CoolingSchedule.Logarithmic,
    "Sawtooth": CoolingSchedule.Sawtooth,
}

DERIVED_MAP: dict[str, tuple[str, str]] = {
    "Chemical potential μ": ("mu", r"chemical potential $\mu / |t|$"),
    "Entropy S": ("entropy", r"entropy $S$"),
    "Heat capacity C": ("heat_capacity", r"heat capacity $C$"),
    "Energy E": ("energy", r"energy $E$"),
    "Free energy F": ("free_energy", r"free energy $F$"),
    "Magnetization M": ("magnetization", r"magnetization $M$"),
    "n_kd": ("nkd", r"$n_{kd}$"),
}


def _solve_mu_from_n(
    *,
    e_k: np.ndarray,
    distro: str,
    n: float,
    kbt: float,
    h_a: float,
    u: float,
) -> tuple[float | None, str | None]:
    """Solve μ for band filling *n* (same Brent defaults as ``ScientificPlotter``)."""
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


def _sidebar_distros(*, key_prefix: str = "") -> list[DistType]:
    st.sidebar.subheader("Distributions")
    names = list(DIST_MAP.keys())
    picked = st.sidebar.multiselect(
        "Select distributions",
        options=names,
        default=["FD", "NFL", "SSL"],
        key=f"{key_prefix}distros",
    )
    return [DIST_MAP[k] for k in picked]


def _sidebar_lattice(*, key_prefix: str = "") -> tuple[str, int, float, float]:
    st.sidebar.subheader("Lattice")
    dim = DIM
    N = int(
        st.sidebar.number_input(
            "N (k-grid per axis)",
            min_value=10,
            max_value=4096,
            value=100,
            step=10,
            key=f"{key_prefix}N",
        )
    )
    t = float(st.sidebar.number_input("t (NN hopping)", value=-1.0, key=f"{key_prefix}t"))
    t_prim = float(st.sidebar.number_input("t' (NNN hopping)", value=0.0, key=f"{key_prefix}t_prim"))
    return dim, N, t, t_prim


def _solver_controls(container: object) -> tuple[str, object | None, float, float, float]:
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


def _df_for_quantity(*, quantity_name: str) -> tuple[pd.DataFrame, str, str]:
    y_col, y_label = DERIVED_MAP[quantity_name]
    if y_col == "mu":
        return st.session_state["mu_df"], y_col, y_label
    if y_col == "entropy":
        return st.session_state["entropy_df"], y_col, y_label
    if y_col == "heat_capacity":
        return st.session_state["heat_cap_df"], y_col, y_label
    if y_col in ("energy", "free_energy"):
        return st.session_state["energy_df"], y_col, y_label
    if y_col == "magnetization":
        return st.session_state["mag_df"], y_col, y_label
    if y_col == "nkd":
        return st.session_state["nkd_df"], y_col, y_label
    raise ValueError(f"Unsupported quantity: {quantity_name}")


def _plot_style_controls() -> dict[str, object]:
    st.subheader("Plot customization")
    c1, c2, c3 = st.columns(3)
    fig_w = float(
        c1.number_input(
            "Figure width (per panel)", min_value=4.0, max_value=24.0, value=10.0, step=0.5
        )
    )
    fig_h = float(
        c2.number_input(
            "Figure height (per panel)", min_value=3.0, max_value=14.0, value=6.0, step=0.5
        )
    )
    cols = int(c3.number_input("Panels per row", min_value=1, max_value=8, value=2, step=1))

    xl_enabled = st.checkbox("Set X axis limits", value=False)
    xlims: tuple[float, float] | None = None
    if xl_enabled:
        x1, x2 = st.columns(2)
        xmin = float(x1.number_input("x min", value=0.0))
        xmax = float(x2.number_input("x max", value=1.0))
        if xmax <= xmin:
            st.warning("x max should be greater than x min.")
        else:
            xlims = (xmin, xmax)

    yl_enabled = st.checkbox("Set Y axis limits", value=False)
    ylims: tuple[float, float] | None = None
    if yl_enabled:
        y1, y2 = st.columns(2)
        ymin = float(y1.number_input("y min", value=0.0))
        ymax = float(y2.number_input("y max", value=1.0))
        if ymax <= ymin:
            st.warning("y max should be greater than y min.")
        else:
            ylims = (ymin, ymax)

    return {
        "figsize": (fig_w, fig_h),
        "cols": cols,
        "xlims": xlims,
        "ylims": ylims,
    }


def _plot_occupancy_map(  # noqa: PLR0913
    *,
    e_k: np.ndarray,
    kx: np.ndarray,
    ky: np.ndarray | None,
    dim: str,
    n: float,
    mu: float,
    kbt: float,
    h_a: float,
    u: float,
    dist_type: DistType,
    distro_label: str,
    spin_choice: str,
) -> None:
    is_nkd = "nkd" in spin_choice
    spin_shift = -1 if "up" in spin_choice else 1
    delta = 0.5 * 2.0 * 1.0 * float(h_a)

    occ = distribution(
        e_k,
        mu=float(mu),
        kBT=float(kbt),
        delta=delta,
        spin_shift=spin_shift,
        dist_type=dist_type,
        U=float(u),
        n_kd=is_nkd,
    )
    plt.close("all")
    c1, c2 = st.columns(2)
    w = float(c1.slider("Chart width", 4.0, 28.0, 10.0, 0.5, key="dist_map_w"))
    h = float(c2.slider("Chart height", 3.0, 20.0, 7.0, 0.5, key="dist_map_h"))

    show_mode = st.radio(
        "Display mode",
        ["2D Map", "Horizontal Cut (k_y=0)", "Diagonal Cut"],
        index=0,
        horizontal=True,
        key="dist_map_mode",
    )

    lim2d = st.checkbox("Custom color/Y limits", value=False, key="dist_map_lim")
    vmin = vmax = None
    if lim2d:
        cl1, cl2 = st.columns(2)
        vmin = float(cl1.number_input("Min", value=0.0, key="dist_map_vmin"))
        vmax = float(cl2.number_input("Max", value=1.0, key="dist_map_vmax"))

    fig, ax = plt.subplots(figsize=(w, h))
    if dim == "2D" and ky is not None and show_mode == "2D Map":
        im = ax.imshow(
            occ,
            extent=[-np.pi, np.pi, -np.pi, np.pi],
            origin="lower",
            aspect="auto",
            vmin=vmin,
            vmax=vmax,
        )
        fig.colorbar(im, ax=ax, label=spin_choice)
        ax.set_xlabel(r"$k_x$")
        ax.set_ylabel(r"$k_y$")
    else:
        if dim == "2D" and ky is not None and show_mode == "Horizontal Cut (k_y=0)":
            mid_idx = occ.shape[0] // 2
            y_vals = occ[mid_idx, :]
            ax.set_xlabel(r"$k_x$ (at $k_y=0$)")
        elif dim == "2D" and ky is not None:
            y_vals = np.diag(occ)
            ax.set_xlabel(r"Diagonal $k$")
        else:
            y_vals = np.diag(occ) if occ.ndim == DIM_2D else occ
            ax.set_xlabel(r"$k$")

        ax.plot(
            kx if kx is not None else np.linspace(-np.pi, np.pi, len(y_vals)),
            y_vals,
            lw=2,
            label=spin_choice,
        )
        ax.set_ylabel("occupancy")
        if vmin is not None and vmax is not None:
            ax.set_ylim(vmin, vmax)
        ax.legend()

    ax.set_title(
        f"{distro_label}, n={n:.6g}, μ={mu:.6g}, kBT={kbt:.6g}, "
        f"H_a={h_a:.6g}, U={u:.6g}"
    )
    st.pyplot(fig, clear_figure=True)


def _view_panels() -> None:  # noqa: PLR0912, PLR0915
    distros = _sidebar_distros(key_prefix="panel_")
    dim, N, t, t_prim = _sidebar_lattice(key_prefix="panel_")

    st.sidebar.subheader("Experiment parameters")
    n_vals, kBT_vals, H_a_vals, U_vals = experiment_parameter_sweeps(st.sidebar)

    with st.sidebar.expander("Solver parameters", expanded=False) as solver_exp:
        solver, solver_opts, a0, b0, x0 = _solver_controls(solver_exp)

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
            e_k, kx, ky = build_dispersion(dim, N, t, t_prim)
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

            st.session_state["mu_df"] = mu_df
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
            with st.spinner("Computing S, C, E, F, M, n_kd…"):
                st.session_state["entropy_df"] = compute_entropy(mu_df, e_k)
                st.session_state["heat_cap_df"] = compute_heat_capacity(mu_df, e_k)
                st.session_state["energy_df"] = compute_energy(mu_df, e_k)
                st.session_state["mag_df"] = calculate_magnetization(mu_df, e_k)
                st.session_state["nkd_df"] = compute_nkd(mu_df, e_k)
            st.success("Finished.")
        except Exception as e:
            prog_ph.empty()
            st.error(str(e))
            st.code(traceback.format_exc())
            return

    mu_df = st.session_state.get("mu_df")
    if mu_df is None:
        st.info("Configure the sidebar and click **Run experiment** to plot sweep panels.")
        return

    if "entropy_df" not in st.session_state:
        st.error("Derived quantities missing — run the experiment again.")
        return

    with st.expander("Plot style", expanded=False):
        style = _plot_style_controls()

    meta = st.session_state["meta"]
    distros = meta["distros"]
    n_arr = meta["n_vals"]
    kBT_arr = meta["kBT_vals"]
    H_arr = meta["H_a_vals"]
    U_arr = meta["U_vals"]

    st.header("Sweep panels")
    q_name = st.selectbox(
        "Quantity", options=list(DERIVED_MAP.keys()), index=0, key="panel_qty"
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

    dfq, y_col, y_label = _df_for_quantity(quantity_name=str(q_name))
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
    var_inline = {"n": r"$n$", "k_BT": r"$k_B T/|t|$", "H_a": r"$H_a$", "U": r"$U/|t|$"}
    fixed_md = (
        ", ".join([f"{var_inline.get(k, k)}={v:.6g}" for k, v in fixed.items()]) or "(none)"
    )
    x_md = var_inline.get(str(x_var), str(x_var))
    pan_md = var_inline.get(str(panel_var_raw), "none") if panel_var_raw != "(none)" else "none"
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


def _view_distributions() -> None:  # noqa: PLR0912, PLR0915
    st.sidebar.caption("Standalone — no experiment run required.")
    dim, N, t, t_prim = _sidebar_lattice(key_prefix="dist_")
    distros = _sidebar_distros(key_prefix="dist_")

    st.sidebar.subheader("State point")
    n_val = float(
        st.sidebar.number_input(
            r"Band filling $n$",
            min_value=0.001,
            max_value=1.0,
            value=0.5,
            step=0.01,
            key="dist_n",
        )
    )
    kbt = float(st.sidebar.number_input(r"$k_B T / |t|$", value=0.1, key="dist_kbt"))
    h_a = float(st.sidebar.number_input(r"$H_a$ (T)", value=0.0, key="dist_ha"))
    u_val = float(st.sidebar.number_input(r"$U / |t|$", value=0.0, key="dist_u"))

    if not distros:
        st.info("Select at least one distribution in the sidebar.")
        return

    e_k, kx, ky = build_dispersion(dim, N, t, t_prim)
    distro_names = [d.name if hasattr(d, "name") else str(d).split(".")[-1] for d in distros]

    st.header("Distributions")
    st.caption(
        "Occupancy and BZ cuts from lattice + $(n, k_B T, H_a, U)$; "
        r"$\mu$ is solved from $n$ for each distribution. "
        "Does not use panel experiment results."
    )

    if dim != "2D":
        st.warning(
            "1D lattice: map/cut view only. For multi-curve 1D BZ cuts with several "
            "$(n, k_B T)$ points, switch to **2D** in the sidebar."
        )

    st.subheader("Occupancy map / 1D cut")
    d_pick = st.selectbox("Distribution", options=distro_names, index=0, key="dist_d_pick")
    mu, mu_err = _solve_mu_from_n(
        e_k=e_k,
        distro=str(d_pick),
        n=n_val,
        kbt=kbt,
        h_a=h_a,
        u=u_val,
    )
    if mu is None:
        st.error(f"Could not solve μ for **{d_pick}**: {mu_err}")
    else:
        st.caption(
            f"Solved μ / |t| = **{mu:.6g}** for n = {n_val:.6g}, "
            f"k_BT = {kbt:.6g}, H_a = {h_a:.6g}, U = {u_val:.6g}."
        )
        spin_choice = st.radio(
            "Spin channel",
            ["up (+1)", "down (-1)", "nkd"],
            index=0,
            horizontal=True,
            key="dist_spin",
        )
        _plot_occupancy_map(
            e_k=e_k,
            kx=kx,
            ky=ky,
            dim=dim,
            n=n_val,
            mu=mu,
            kbt=kbt,
            h_a=h_a,
            u=u_val,
            dist_type=DIST_MAP[str(d_pick)],
            distro_label=str(d_pick),
            spin_choice=spin_choice,
        )

    if dim == "2D" and ky is not None:
        st.divider()
        st.subheader("1D BZ cuts")
        st.caption(
            "Uses `ScientificPlotter` — each curve fixes **n** and **$k_B T$** "
            "(μ solved per point). Independent of the Panels experiment."
        )
        slice_type = st.selectbox("1D cut", ["horizontal", "diagonal"], index=0, key="dist_slice")
        n_text = st.text_input(
            "n values (comma-separated)",
            value="0.1, 0.5, 0.9",
            key="dist_n_list",
        )
        kbt_text = st.text_input(
            "$k_BT$ values (comma-separated)",
            value="0.01, 0.1",
            key="dist_kbt_list",
        )
        n_pick = parse_float_list(n_text) or [0.1, 0.5, 0.9]
        k_pick = parse_float_list(kbt_text) or [0.01, 0.1]
        H_plot = float(st.number_input("$H_a$ for 1D plot", value=h_a, key="dist_h_plot"))
        U_plot = float(st.number_input("$U$ for 1D plot", value=u_val, key="dist_u_plot"))
        n_kd_flag = st.checkbox(
            "n_kd channel / correlator view (Hubbard tooling)", value=False, key="dist_nkd"
        )
        y0, y1 = st.slider("y axis (occupancy)", -0.05, 2.2, (-0.05, 2.2), key="dist_ylim")
        fig_w = st.slider("Figure width (1D)", 10, 48, 30, key="dist_fig_w")
        fig_h = st.slider("Figure height (1D)", 2, 8, 4, key="dist_fig_h")

        if not n_pick or not k_pick:
            st.error("Enter at least one **n** and one **$k_BT$** value.")
        else:
            plt.close("all")
            plotter = ScientificPlotter()
            try:
                plotter.plot_energy_distribution_1d(
                    e_k=e_k,
                    slice_type=slice_type,
                    dim=2,
                    kx=kx,
                    n=[float(x) for x in n_pick],
                    kBTs=[float(x) for x in k_pick],
                    H=H_plot,
                    U=float(U_plot),
                    n_kd=n_kd_flag,
                    distributions=distro_names,
                    ylims=(float(y0), float(y1)),
                    figsize=(float(fig_w), float(fig_h)),
                )
                show_current_figure(st)
                st.caption(
                    f"Distributions: n={len(n_pick)}, kBT={len(k_pick)}, "
                    f"H_a={H_plot:.6g}, U={U_plot:.6g}, slice={slice_type}."
                )
            except Exception as e:
                st.error(f"Error drawing 1D cuts: {e}")


def main() -> None:
    st.title("Quantum statistics solver")
    st.caption(
        "Notebook-equivalent sweeps and BZ distributions — "
        f"project root: `{SOLVER_ROOT}`"
    )

    view = st.sidebar.radio(
        "View",
        ["Panels", "Distributions"],
        index=0,
        horizontal=True,
        key="main_view",
    )

    if view == "Panels":
        _view_panels()
    else:
        _view_distributions()


if __name__ == "__main__":
    main()
