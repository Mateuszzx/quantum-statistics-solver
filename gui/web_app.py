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
from gui.panel_plot import mu_for_point, panel_plot
from gui.paths import quantum_statistics_solver_root
from gui.scientific_plot import ScientificPlotter
from gui.solver_options import build_brent_options, build_sa_options
from gui.streamlit_utils import show_current_figure, st_progress_pair
from gui.sweep_grids import build_1d_sweep
from qss.distributions import DistType, distribution
from qss.energy import energy as compute_energy
from qss.entropy import entropy as compute_entropy
from qss.magnetization import calculate_magnetization
from qss.nkd import n_kd as compute_nkd
from qss.solver_options import CoolingSchedule
from utils.experiment_utils import run_experiment

st.set_page_config(page_title="Quantum Statistics Solver", layout="wide")

SOLVER_ROOT = quantum_statistics_solver_root()

DIM_2D = 2
MIN_SWEEP_COUNT = 2

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
    "Energy E": ("energy", r"energy $E$"),
    "Free energy F": ("free_energy", r"free energy $F$"),
    "Magnetization M": ("magnetization", r"magnetization $M$"),
    "n_kd": ("nkd", r"$n_{kd}$"),
}


def _sidebar_distros() -> list[DistType]:
    st.sidebar.subheader("Distributions")
    names = list(DIST_MAP.keys())
    picked = st.sidebar.multiselect(
        "Select distributions", options=names, default=["FD", "NFL", "SSL"]
    )
    return [DIST_MAP[k] for k in picked]


def _sidebar_lattice() -> tuple[str, int, float, float]:
    st.sidebar.subheader("Lattice")
    dim = st.sidebar.radio("Dimensionality", ["2D", "1D"], index=0, horizontal=True)
    N = int(
        st.sidebar.number_input(
            "N (k-grid per axis)", min_value=10, max_value=4096, value=100, step=10
        )
    )
    t = float(st.sidebar.number_input("t (NN hopping)", value=-1.0))
    t_prim = float(st.sidebar.number_input("t' (NNN hopping)", value=0.0))
    return dim, N, t, t_prim


def _sidebar_sweeps() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    st.sidebar.subheader("Sweeps")
    c1, c2 = st.sidebar.columns(2)
    n_min = c1.number_input("n min", 0.001, 1.0, 0.1)
    n_max = c2.number_input("n max", 0.001, 1.0, 0.99)
    n_count = int(st.sidebar.number_input("n points", 2, 500, 20))

    st.sidebar.markdown(r"$$k_B T/|t|$$")
    k_mode = st.sidebar.radio("$k_BT$ spacing", ["linear", "log"], horizontal=True, key="kbt_mode")
    k1, k2 = st.sidebar.columns(2)
    k_min = float(k1.number_input("$k_BT$ min", value=0.01, format="%.1e", key="kmin"))
    k_max = float(k2.number_input("$k_BT$ max", value=1.0, format="%f", key="kmax"))
    k_count = int(st.sidebar.number_input("$k_BT$ points", 1, 200, 4, key="kcn"))
    if k_mode == "log" and (k_min <= 0 or k_max <= 0):
        st.sidebar.error("Log $k_BT$ sweep needs positive min and max.")
    kBT_vals = build_1d_sweep(
        k_min, k_max, max(2, k_count), "linear" if k_mode == "linear" else "log"
    )

    st.sidebar.markdown(r"$$H_a\ \mathrm{(T)}$$")
    h_mode = st.sidebar.radio("$H_a$ spacing", ["linear", "log"], horizontal=True, key="hm")
    h1, h2 = st.sidebar.columns(2)
    h_min = float(h1.number_input("$H_a$ min", value=0.0, format="%.1e", key="hmin"))
    h_max = float(h2.number_input("$H_a$ max", value=0.0, key="hmax"))
    h_count = int(st.sidebar.number_input("$H_a$ points", 1, 200, 1, key="hcn"))
    if h_mode == "log" and h_count >= MIN_SWEEP_COUNT and (h_min <= 0 or h_max <= 0):
        st.sidebar.error("Log H_a sweep needs strictly positive min and max.")
    H_a_vals = build_1d_sweep(h_min, h_max, h_count, "linear" if h_mode == "linear" else "log")

    n_vals = build_1d_sweep(float(n_min), float(n_max), n_count, "linear")
    return n_vals, kBT_vals, H_a_vals


def _sidebar_physics() -> float:
    st.sidebar.subheader("Hubbard")
    st.sidebar.markdown(r"$$U/|t|$$")
    return float(st.sidebar.number_input("U / |t|", value=0.0))


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
        with container.expander("Simulated annealing", expanded=False):
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
        with container.expander("Brent / bracketing", expanded=False):
            abs_tol = container.number_input("abs_tol", value=1e-9, format="%.1e", key="br_abs")
            rel_tol = container.number_input("rel_tol", value=1e-6, format="%.1e", key="br_rel")
            max_iters = container.number_input("max_iters", value=5000, key="br_max_iters")
            expand_factor = container.number_input("expand_factor", value=2.0, key="br_expand")
            max_expands = container.number_input("max_expands", value=164, key="br_max_expands")
            initial_step = container.number_input("initial_step", value=1.0, key="br_initial")
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


def main() -> None:  # noqa: PLR0912, PLR0915
    st.title("Quantum statistics solver")
    st.caption(
        "Notebook-equivalent sweeps over band filling, temperature, and field — "
        f"project root: `{SOLVER_ROOT}`"
    )

    distros = _sidebar_distros()
    dim, N, t, t_prim = _sidebar_lattice()
    n_vals, kBT_vals, H_a_vals = _sidebar_sweeps()
    U = _sidebar_physics()

    with st.expander("Solver parameters", expanded=False):
        solver, solver_opts, a0, b0, x0 = _solver_controls(st)

    logs_dir = SOLVER_ROOT / "logs"
    out_csv = SOLVER_ROOT / "data" / "gui_results.csv"
    logs_dir.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    est = len(distros) * len(n_vals) * len(kBT_vals) * len(H_a_vals)
    st.sidebar.metric("Total solves (approx.)", f"{est:,}")

    if st.sidebar.button("Run experiment", type="primary"):
        if not distros:
            st.error("Select at least one distribution.")
            return
        prog_cb, prog_ph = st_progress_pair(st)
        try:
            e_k, kx, ky = build_dispersion(dim, N, t, t_prim)
            n_run = n_vals
            kbt_run = kBT_vals
            h_run = H_a_vals
            out_csv_run = out_csv
            with st.spinner("Running chemical potential solves…"):
                mu_df = run_experiment(
                    n_vals=n_run,
                    kBT_vals=kbt_run,
                    h_a_vals=h_run,
                    distros=distros,
                    e_k=e_k,
                    N=int(e_k.shape[0]),
                    solver=solver,
                    solver_opts=solver_opts,
                    logs_dir=str(logs_dir),
                    output_file=str(out_csv_run),
                    x0=x0,
                    a0=a0,
                    b0=b0,
                    U=U,
                    progress=prog_cb,
                )
            prog_ph.empty()
            if mu_df.empty:
                st.error("No results returned. Check stderr / logs for solver errors.")
                return

            st.session_state["mu_df"] = mu_df
            st.session_state["e_k"] = e_k
            st.session_state["kx"] = kx
            st.session_state["ky"] = ky
            st.session_state["meta"] = {
                "n_vals": np.array(n_run, dtype=float),
                "kBT_vals": np.array(kbt_run, dtype=float),
                "H_a_vals": np.array(h_run, dtype=float),
                "distros": distros,
                "t_prim": t_prim,
                "N": int(e_k.shape[0]),
                "dim": dim,
                "U": U,
            }
            with st.spinner("Computing S, E, F, M, n_kd…"):
                st.session_state["entropy_df"] = compute_entropy(mu_df, e_k, U=U)
                st.session_state["energy_df"] = compute_energy(mu_df, e_k, U=U)
                st.session_state["mag_df"] = calculate_magnetization(mu_df, e_k, U=U)
                st.session_state["nkd_df"] = compute_nkd(mu_df, e_k, U=U)
            st.success("Finished.")
        except Exception as e:
            prog_ph.empty()
            st.error(str(e))
            st.code(traceback.format_exc())
            return

    mu_df = st.session_state.get("mu_df")
    if mu_df is None:
        st.info("Configure the sidebar and click **Run experiment**.")
        return

    if mu_df is not None:
        meta = st.session_state["meta"]
        e_k = st.session_state["e_k"]
        Umeta = float(meta["U"])
        if "entropy_df" not in st.session_state:
            with st.spinner("Computing S, E, F, M, n_kd…"):
                st.session_state["entropy_df"] = compute_entropy(mu_df, e_k, U=Umeta)
                st.session_state["energy_df"] = compute_energy(mu_df, e_k, U=Umeta)
                st.session_state["mag_df"] = calculate_magnetization(mu_df, e_k, U=Umeta)
                st.session_state["nkd_df"] = compute_nkd(mu_df, e_k, U=Umeta)

    with st.expander("Plot style", expanded=False):
        style = _plot_style_controls()

    meta = st.session_state["meta"]
    kx = st.session_state["kx"]
    ky = st.session_state["ky"]
    distros = meta["distros"]
    n_arr = meta["n_vals"]
    kBT_arr = meta["kBT_vals"]
    H_arr = meta["H_a_vals"]

    tab_panels, tab_dist = st.tabs(["Panels", "Distributions"])

    with tab_panels, st.expander("Panels", expanded=True):
        q_name = st.selectbox(
            "Quantity", options=list(DERIVED_MAP.keys()), index=0, key="panel_qty"
        )
        x_var = st.selectbox(
            "X axis variable", options=["n", "k_BT", "H_a"], index=0, key="panel_x"
        )
        x_scale = st.radio(
            "X scale", ["linear", "log"], horizontal=True, index=0, key="panel_xscale"
        )
        panel_var_raw = st.selectbox(
            "Sweep in panels", options=["(none)", "k_BT", "H_a"], index=1, key="panel_pan"
        )
        panel_var = None if panel_var_raw == "(none)" else str(panel_var_raw)
        if panel_var == x_var:
            st.warning("Panel variable must be different from X axis. Switching panels to (none).")
            panel_var = None

        n_cols = int(style["cols"])

        fixed: dict[str, float] = {}
        for var, opts in (("n", n_arr), ("k_BT", kBT_arr), ("H_a", H_arr)):
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
        var_inline = {"n": r"$n$", "k_BT": r"$k_B T/|t|$", "H_a": r"$H_a$"}
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

    with tab_dist:
        st.subheader("1D Distribution")
        if meta["dim"] != "2D":
            st.warning(
                "Advanced 1D cuts are designed for a 2D lattice. Please select 2D in the sidebar."
            )
        else:
            slice_type = st.selectbox("1D cut", ["horizontal", "diagonal"], index=0)
            n_opts = list(n_arr.astype(float))
            k_opts = list(kBT_arr.astype(float))
            n_pick = st.multiselect(
                "n values", options=n_opts, default=n_opts[: min(3, len(n_opts))]
            )
            k_pick = st.multiselect(
                "$k_BT$ values", options=k_opts, default=k_opts[: min(2, len(k_opts))]
            )
            if not n_pick or not k_pick:
                st.error("Select at least one value for **n** and for **$k_BT$**.")
            else:
                H_plot = float(st.selectbox("H_a for plot", list(H_arr.astype(float))))
                n_kd_flag = st.checkbox(
                    "n_kd channel / correlator view (Hubbard tooling)", value=False
                )
                y0, y1 = st.slider("y axis (occupancy)", -0.05, 2.2, (-0.05, 2.2))
                fig_w = st.slider("Figure width (1D)", 10, 48, 30)
                fig_h = st.slider("Figure height (1D)", 2, 8, 4)

                names = [d.name if hasattr(d, "name") else str(d) for d in distros]
                plt.close("all")

                plotter = ScientificPlotter()
                try:
                    plotter.plot_energy_distribution_1d(
                        e_k=e_k,
                        slice_type=slice_type,
                        dim=2 if meta["dim"] == "2D" else 1,
                        kx=kx,
                        n=[float(x) for x in n_pick],
                        kBTs=[float(x) for x in k_pick],
                        H=H_plot,
                        U=float(meta["U"]),
                        n_kd=n_kd_flag,
                        distributions=names,
                        ylims=(float(y0), float(y1)),
                        figsize=(float(fig_w), float(fig_h)),
                    )
                    show_current_figure(st)
                    st.caption(
                        f"Dystrybucje: n={len(n_pick)} punktow, kBT={len(k_pick)} punktow, "
                        f"H_a={H_plot:.6g}, U={float(meta['U']):.6g}, slice={slice_type}."
                    )
                except Exception as e:
                    st.error(f"Error drawing 1D cuts: {e}")

        if meta["dim"] == "2D" and ky is not None:
            st.divider()
            st.subheader("2D occupancy map (uses solved μ)")
            mu_df = st.session_state["mu_df"]
            distro_names = [
                d.name if hasattr(d, "name") else str(d).split(".")[-1] for d in distros
            ]
            d_pick = st.selectbox("Distribution", options=distro_names, index=0, key="occ2d_d")
            n2 = float(st.selectbox("n", options=list(n_arr.astype(float)), key="occ2d_n"))
            k2 = float(st.selectbox("$k_BT$", options=list(kBT_arr.astype(float)), key="occ2d_k"))
            h2 = float(st.selectbox("H_a", options=list(H_arr.astype(float)), key="occ2d_h"))
            mu_val = mu_for_point(mu_df, distro=str(d_pick), n=n2, kBT=k2, H_a=h2)
            if mu_val is None:
                st.error("Could not find μ for this point in results.")
            else:
                st.caption(f"Using μ = {mu_val:.6g}")
                delta = 0.5 * 2.0 * 1.0 * float(h2)
                occ = distribution(
                    e_k,
                    mu=float(mu_val),
                    kBT=float(k2),
                    delta=delta,
                    spin_shift=1,
                    dist_type=DIST_MAP[str(d_pick)],
                    U=float(meta["U"]),
                )
                plt.close("all")
                c1, c2 = st.columns(2)
                w = float(c1.slider("2D width", 4.0, 28.0, 10.0, 0.5))
                h = float(c2.slider("2D height", 3.0, 20.0, 7.0, 0.5))
                lim2d = st.checkbox("Custom color limits", value=False, key="occ2d_lim")
                vmin = vmax = None
                if lim2d:
                    cl1, cl2 = st.columns(2)
                    vmin = float(cl1.number_input("Color min", value=0.0, key="occ2d_vmin"))
                    vmax = float(cl2.number_input("Color max", value=1.0, key="occ2d_vmax"))
                fig, ax = plt.subplots(figsize=(w, h))
                im = ax.imshow(
                    occ,
                    extent=[-np.pi, np.pi, -np.pi, np.pi],
                    origin="lower",
                    aspect="auto",
                    vmin=vmin,
                    vmax=vmax,
                )
                fig.colorbar(im, ax=ax, label="occupancy (one spin)")
                ax.set_xlabel(r"$k_x$")
                ax.set_ylabel(r"$k_y$")
                ax.set_title(f"{d_pick}, n={n2:.6g}, kBT={k2:.6g}, H_a={h2:.6g}")
                st.pyplot(fig, clear_figure=True)


if __name__ == "__main__":
    main()
