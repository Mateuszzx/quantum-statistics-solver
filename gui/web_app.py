"""
Streamlit GUI for quantum-statistics-solver: reproduces notebook workflows (sweeps in
``n``, ``k_B T``, ``H_a``; μ, entropy, energy, free energy, magnetization, ``n_{kd}``;
distribution / BZ cuts when ``qss-results/tools`` is available).

Run from the ``quantum-statistics-solver`` directory::

    streamlit run gui/web_app.py
"""

from __future__ import annotations

import traceback
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from qss.distributions import DistType, distribution
from qss.energy import energy as compute_energy
from qss.entropy import entropy as compute_entropy
from qss.magnetization import calculate_magnetization
from qss.nkd import n_kd as compute_nkd
from qss.solver_options import CoolingSchedule

import utils.scientific_plot as sp
from utils.experiment_utils import run_experiment
from utils.gui_lattice import build_dispersion
from utils.gui_notebook_plots import load_scientific_plotter_v2, plot_distributions_notebook_style
from utils.gui_paths import quantum_statistics_solver_root
from utils.gui_solver_options import build_brent_options, build_sa_options
from utils.gui_streamlit import show_current_figure, st_progress_pair
from utils.gui_sweep_grids import betas_for_scientific_plot_n, build_1d_sweep

st.set_page_config(page_title="Quantum Statistics Solver", layout="wide")

SOLVER_ROOT = quantum_statistics_solver_root()

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


def _sidebar_distros() -> list[DistType]:
    st.sidebar.subheader("Distributions")
    names = list(DIST_MAP.keys())
    picked = st.sidebar.multiselect("Select distributions", options=names, default=["FD", "NFL", "SSL"])
    return [DIST_MAP[k] for k in picked]


def _sidebar_lattice() -> tuple[str, int, float, float]:
    st.sidebar.subheader("Lattice")
    dim = st.sidebar.radio("Dimensionality", ["2D", "1D"], index=0, horizontal=True)
    N = int(st.sidebar.number_input("N (k-grid per axis)", min_value=10, max_value=4096, value=100, step=10))
    t = float(st.sidebar.number_input("t (NN hopping)", value=-1.0))
    t_prim = float(st.sidebar.number_input("t' (NNN hopping)", value=0.0))
    return dim, N, t, t_prim


def _sidebar_sweeps() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    st.sidebar.subheader("Sweeps (notebook-style)")
    c1, c2 = st.sidebar.columns(2)
    n_min = c1.number_input("n min", 0.001, 1.0, 0.1)
    n_max = c2.number_input("n max", 0.001, 1.0, 0.99)
    n_count = int(st.sidebar.number_input("n points", 2, 500, 20))

    st.sidebar.markdown("**k_B T / |t|**")
    k_mode = st.sidebar.radio("k_BT spacing", ["linear", "log"], horizontal=True, key="kbt_mode")
    k1, k2 = st.sidebar.columns(2)
    k_min = float(k1.number_input("k_BT min", value=0.01, format="%f", key="kmin"))
    k_max = float(k2.number_input("k_BT max", value=1.0, format="%f", key="kmax"))
    k_count = int(st.sidebar.number_input("k_BT points", 1, 200, 5, key="kcn"))
    if k_mode == "log" and (k_min <= 0 or k_max <= 0):
        st.sidebar.error("Log k_BT sweep needs positive min and max.")
    kBT_vals = build_1d_sweep(k_min, k_max, max(2, k_count), "linear" if k_mode == "linear" else "log")

    st.sidebar.markdown("**H_a (T)**")
    h_mode = st.sidebar.radio("H_a spacing", ["linear", "log"], horizontal=True, key="hm")
    h1, h2 = st.sidebar.columns(2)
    h_min = float(h1.number_input("H_a min", value=0.0, key="hmin"))
    h_max = float(h2.number_input("H_a max", value=0.0, key="hmax"))
    h_count = int(st.sidebar.number_input("H_a points", 1, 200, 1, key="hcn"))
    if h_mode == "log" and h_count >= 2 and (h_min <= 0 or h_max <= 0):
        st.sidebar.error("Log H_a sweep needs strictly positive min and max.")
    H_a_vals = build_1d_sweep(h_min, h_max, h_count, "linear" if h_mode == "linear" else "log")

    n_vals = build_1d_sweep(float(n_min), float(n_max), n_count, "linear")
    return n_vals, kBT_vals, H_a_vals


def _sidebar_physics() -> float:
    st.sidebar.subheader("Hubbard")
    return float(st.sidebar.number_input("U / |t|", value=0.0))


def _sidebar_solver() -> tuple[str, object | None, float, float, float]:
    st.sidebar.subheader("Solver")
    solver = st.sidebar.selectbox("Algorithm", ["Brent", "SA"], index=0)
    opts: object | None = None

    st.sidebar.subheader("Search / guess")
    search_min = float(st.sidebar.number_input("Bracket a0", value=-10.0, key="a0"))
    search_max = float(st.sidebar.number_input("Bracket b0", value=10.0, key="b0"))
    use_random = st.sidebar.checkbox("Random initial guess", value=False)
    if use_random:
        x0 = float(np.random.uniform(search_min, search_max))
        st.sidebar.caption(f"Using x0 = {x0:.6g}")
    else:
        x0 = float(st.sidebar.number_input("Initial guess x0", value=0.0))

    if solver == "SA":
        with st.sidebar.expander("Simulated annealing", expanded=False):
            max_iters = st.number_input("Max iterations", value=50000)
            min_iters = st.number_input("Min iterations", value=2000)
            stall_iters = st.number_input("Stall iterations", value=1000)
            T0 = st.number_input("T0", value=1000.0)
            T_min = st.number_input("T_min", value=1e-9, format="%.1e")
            alpha = st.number_input("alpha", value=0.999, min_value=1e-6, max_value=1.0)
            step_decay = st.number_input("step_decay", value=0.999, min_value=1e-6, max_value=1.0)
            step_sigma = st.number_input("step_sigma", value=5.0)
            sn = st.selectbox("Cooling schedule", list(SCHED_MAP.keys()))
            scale_sigma = st.checkbox("Scale sigma by T", value=True)
            reflect_bounds = st.checkbox("Reflect bounds", value=True)
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
        with st.sidebar.expander("Brent / bracketing", expanded=False):
            abs_tol = st.number_input("abs_tol", value=1e-9, format="%.1e")
            rel_tol = st.number_input("rel_tol", value=1e-6, format="%.1e")
            max_iters = st.number_input("max_iters", value=5000)
            expand_factor = st.number_input("expand_factor", value=2.0)
            max_expands = st.number_input("max_expands", value=164)
            initial_step = st.number_input("initial_step", value=1.0)
            brk = st.selectbox("Bracketer", ["Expansion", "Bisection"])
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


def _plot_scalar_vs_n(
    df: pd.DataFrame,
    *,
    n_vals: np.ndarray,
    kBT_vals: np.ndarray,
    distros: list[DistType],
    H_pick: float,
    y_col: str,
    y_label: str,
    t_prim: float,
    N: int,
    show_log: bool = False,
    figsize: tuple[float, float] = (10.0, 6.0),
    xlims: tuple[float, float] | None = None,
    ylims: tuple[float, float] | None = None,
    text_on_plot: bool = True,
    cols: int | None = None,
) -> None:
    betas = betas_for_scientific_plot_n(kBT_vals)
    plt.close("all")
    sp.plot_result_n(
        results=df,
        n=n_vals,
        betas=betas,
        y_label=y_label,
        distros=distros,
        H_a=float(H_pick),
        y_col=y_col,
        legend_position="center left",
        figsize=figsize,
        cols=(min(5, max(1, len(kBT_vals))) if cols is None else int(cols)),
        t_prim=float(t_prim),
        N=int(N),
        show_log_vals=show_log,
        xlims=xlims,
        ylims=ylims,
        text_on_plot=text_on_plot,
    )


def _plot_style_controls() -> dict[str, object]:
    st.subheader("Plot customization")
    c1, c2, c3 = st.columns(3)
    fig_w = float(c1.number_input("Figure width (per panel)", min_value=4.0, max_value=24.0, value=10.0, step=0.5))
    fig_h = float(c2.number_input("Figure height (per panel)", min_value=3.0, max_value=14.0, value=6.0, step=0.5))
    cols = int(c3.number_input("Panels per row", min_value=1, max_value=8, value=5, step=1))

    text_placement = st.radio("Panel text placement", ["On plot", "Below plot"], horizontal=True, index=1)

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
        "text_on_plot": text_placement == "On plot",
    }


def _render_panel_note(style: dict[str, object], *, n_count: int, kbt_count: int, fixed: str) -> None:
    if bool(style["text_on_plot"]):
        return
    st.caption(
        f"t = -1, N panels: {kbt_count}, n points: {n_count}, fixed {fixed}. "
        "Plot metadata is displayed here instead of inside each subplot."
    )


def _extra_sweeps_for_scalar(
    *,
    key_prefix: str,
    y_col: str,
    y_label: str,
    df: pd.DataFrame,
    n_arr: np.ndarray,
    kBT_arr: np.ndarray,
    H_pick: float,
    beta_for_hplot: float,
    distros: list[DistType],
    meta: dict[str, object],
    style: dict[str, object],
    kBT_pick: float,
) -> None:
    if st.checkbox(f"Also plot {y_col} vs k_B T (fixed n, H)", value=False, key=f"{key_prefix}_vs_kbt"):
        n_fix = float(st.selectbox("n (fixed)", list(n_arr.astype(float)), key=f"nf_{key_prefix}"))
        plt.close("all")
        sp.plot_results_temperature(
            df,
            n=np.array([n_fix], dtype=float),
            betas=betas_for_scientific_plot_n(kBT_arr),
            y_label=y_label,
            distros=distros,
            H_a=float(H_pick),
            y_col=y_col,
            cols=int(style["cols"]),
            figsize=style["figsize"],
            t_prim=float(meta["t_prim"]),
            N=int(meta["N"]),
            xlims=style["xlims"],
            ylims=style["ylims"],
            text_on_plot=bool(style["text_on_plot"]),
        )
        show_current_figure(st)
        _render_panel_note(style, n_count=len(n_arr), kbt_count=len(kBT_arr), fixed=f"H_a={H_pick:.6g}")

    if st.checkbox(f"Also plot {y_col} vs H_a (fixed n, k_BT)", value=False, key=f"{key_prefix}_vs_h"):
        n_fix2 = float(st.selectbox("n (fixed) ", list(n_arr.astype(float)), key=f"nf2_{key_prefix}"))
        plt.close("all")
        sp.plot_results_Ha(
            df,
            n=np.array([n_fix2], dtype=float),
            beta=beta_for_hplot,
            y_label=y_label,
            distros=distros,
            y_col=y_col,
            cols=int(style["cols"]),
            figsize=style["figsize"],
            t_prim=float(meta["t_prim"]),
            N=int(meta["N"]),
            xlims=style["xlims"],
            ylims=style["ylims"],
            text_on_plot=bool(style["text_on_plot"]),
        )
        show_current_figure(st)
        _render_panel_note(style, n_count=len(n_arr), kbt_count=len(kBT_arr), fixed=f"kBT={kBT_pick:.6g}")


def main() -> None:
    st.title("Quantum statistics solver")
    st.caption(
        "Notebook-equivalent sweeps over band filling, temperature, and field — "
        f"project root: `{SOLVER_ROOT}`"
    )

    distros = _sidebar_distros()
    dim, N, t, t_prim = _sidebar_lattice()
    n_vals, kBT_vals, H_a_vals = _sidebar_sweeps()
    U = _sidebar_physics()
    solver, solver_opts, a0, b0, x0 = _sidebar_solver()

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
                "n_vals": n_vals,
                "kBT_vals": kBT_vals,
                "H_a_vals": H_a_vals,
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
        except Exception as e:  # noqa: BLE001
            prog_ph.empty()
            st.error(str(e))
            st.code(traceback.format_exc())
            return

    mu_df = st.session_state.get("mu_df")
    if mu_df is None:
        st.info("Configure the sidebar and click **Run experiment**.")
        return

    meta = st.session_state["meta"]
    e_k = st.session_state["e_k"]
    Umeta = float(meta["U"])
    if "entropy_df" not in st.session_state:
        with st.spinner("Computing S, E, F, M, n_kd…"):
            st.session_state["entropy_df"] = compute_entropy(mu_df, e_k, U=Umeta)
            st.session_state["energy_df"] = compute_energy(mu_df, e_k, U=Umeta)
            st.session_state["mag_df"] = calculate_magnetization(mu_df, e_k, U=Umeta)
            st.session_state["nkd_df"] = compute_nkd(mu_df, e_k, U=Umeta)

    kx = st.session_state["kx"]
    ky = st.session_state["ky"]
    distros = meta["distros"]
    n_arr = meta["n_vals"]
    kBT_arr = meta["kBT_vals"]
    H_arr = meta["H_a_vals"]

    H_pick = st.selectbox(
        "Fixed H_a for “vs n” panels (matches ``plot_result_n``)",
        options=list(H_arr.astype(float)),
        format_func=lambda x: f"{x:.6g}",
    )
    kBT_pick = st.selectbox(
        "Fixed k_B T / |t| for “vs H_a” plots",
        options=list(kBT_arr.astype(float)),
        format_func=lambda x: f"{x:.6g}",
    )
    beta_for_hplot = float(1.0 / float(kBT_pick))
    with st.expander("Plot style", expanded=False):
        style = _plot_style_controls()

    tab_mu, tab_s, tab_e, tab_f, tab_m, tab_nkd, tab_dist = st.tabs(
        ["μ", "Entropy S", "Energy E", "Free energy F", "Magnetization M", "n_kd", "Distributions"]
    )

    with tab_mu:
        _plot_scalar_vs_n(
            mu_df,
            n_vals=n_arr,
            kBT_vals=kBT_arr,
            distros=distros,
            H_pick=H_pick,
            y_col="mu",
            y_label=r"chemical potential $\mu / |t|$",
            t_prim=meta["t_prim"],
            N=meta["N"],
            figsize=style["figsize"],
            xlims=style["xlims"],
            ylims=style["ylims"],
            text_on_plot=bool(style["text_on_plot"]),
            cols=int(style["cols"]),
        )
        show_current_figure(st)
        _render_panel_note(style, n_count=len(n_arr), kbt_count=len(kBT_arr), fixed=f"H_a={H_pick:.6g}")
        with st.expander("Table"):
            st.dataframe(mu_df.reset_index())
            st.download_button(
                "Download μ CSV",
                mu_df.reset_index().to_csv(index=False).encode("utf-8"),
                "gui_mu.csv",
                "text/csv",
            )
        if st.checkbox("Also plot μ vs k_B T (fixed n, H)", value=False, key="mu_vs_kbt"):
            n_fix = float(st.selectbox("n (fixed)", list(n_arr.astype(float)), key="nf_mu"))
            plt.close("all")
            sp.plot_results_temperature(
                mu_df,
                n=np.array([n_fix], dtype=float),
                betas=betas_for_scientific_plot_n(kBT_arr),
                y_label=r"$\mu / |t|$",
                distros=distros,
                H_a=float(H_pick),
                y_col="mu",
                cols=int(style["cols"]),
                figsize=style["figsize"],
                t_prim=float(meta["t_prim"]),
                N=meta["N"],
                xlims=style["xlims"],
                ylims=style["ylims"],
                text_on_plot=bool(style["text_on_plot"]),
            )
            show_current_figure(st)
            _render_panel_note(style, n_count=len(n_arr), kbt_count=len(kBT_arr), fixed=f"H_a={H_pick:.6g}")
        if st.checkbox("Also plot μ vs H_a (fixed n, k_BT)", value=False, key="mu_vs_h"):
            n_fix2 = float(st.selectbox("n (fixed) ", list(n_arr.astype(float)), key="nf2_mu"))
            plt.close("all")
            sp.plot_results_Ha(
                mu_df,
                n=np.array([n_fix2], dtype=float),
                beta=beta_for_hplot,
                y_label=r"$\mu / |t|$",
                distros=distros,
                y_col="mu",
                cols=int(style["cols"]),
                figsize=style["figsize"],
                t_prim=float(meta["t_prim"]),
                N=meta["N"],
                xlims=style["xlims"],
                ylims=style["ylims"],
                text_on_plot=bool(style["text_on_plot"]),
            )
            show_current_figure(st)
            _render_panel_note(style, n_count=len(n_arr), kbt_count=len(kBT_arr), fixed=f"kBT={kBT_pick:.6g}")

    def _derived(
        df: pd.DataFrame,
        y_col: str,
        y_label: str,
        *,
        show_log: bool = False,
        with_extra_sweeps: bool = False,
        key_prefix: str = "derived",
    ) -> None:
        _plot_scalar_vs_n(
            df,
            n_vals=n_arr,
            kBT_vals=kBT_arr,
            distros=distros,
            H_pick=H_pick,
            y_col=y_col,
            y_label=y_label,
            t_prim=meta["t_prim"],
            N=meta["N"],
            show_log=show_log,
            figsize=style["figsize"],
            xlims=style["xlims"],
            ylims=style["ylims"],
            text_on_plot=bool(style["text_on_plot"]),
            cols=int(style["cols"]),
        )
        show_current_figure(st)
        _render_panel_note(style, n_count=len(n_arr), kbt_count=len(kBT_arr), fixed=f"H_a={H_pick:.6g}")
        if with_extra_sweeps:
            _extra_sweeps_for_scalar(
                key_prefix=key_prefix,
                y_col=y_col,
                y_label=y_label,
                df=df,
                n_arr=n_arr,
                kBT_arr=kBT_arr,
                H_pick=H_pick,
                beta_for_hplot=beta_for_hplot,
                distros=distros,
                meta=meta,
                style=style,
                kBT_pick=float(kBT_pick),
            )
        with st.expander("Table"):
            st.dataframe(df.reset_index())
            st.download_button(
                f"Download {y_col} CSV",
                df.reset_index().to_csv(index=False).encode("utf-8"),
                f"gui_{y_col}.csv",
                "text/csv",
            )

    with tab_s:
        _derived(
            st.session_state["entropy_df"],
            "entropy",
            r"entropy $S$",
            show_log=True,
            with_extra_sweeps=True,
            key_prefix="entropy",
        )

    with tab_e:
        _derived(
            st.session_state["energy_df"],
            "energy",
            r"energy $E$",
            with_extra_sweeps=True,
            key_prefix="energy",
        )

    with tab_f:
        _derived(st.session_state["energy_df"], "free_energy", r"free energy $F$")

    with tab_m:
        _derived(st.session_state["mag_df"], "magnetization", r"magnetization $M$")

    with tab_nkd:
        _derived(st.session_state["nkd_df"], "nkd", r"$n_{kd}$")

    with tab_dist:
        st.subheader("BZ occupancy (notebook-style)")
        PlotterCls = load_scientific_plotter_v2()
        if PlotterCls is None:
            st.warning(
                "Advanced cuts need the parent **qss-results** checkout with `tools/scientific_plot_v2.py` "
                "next to this submodule (same layout as the Jupyter notebooks)."
            )
            st.markdown("**Quick view** (single ``μ`` guess — qualitative only).")
            mu_guess = st.number_input("μ guess", value=0.0)
            kbt0 = float(st.selectbox("k_BT", list(kBT_arr.astype(float))))
            H0 = float(st.selectbox("H_a", list(H_arr.astype(float))))
            d0 = distros[0]
            delta = 0.5 * 2.0 * 1.0 * H0
            occ = distribution(
                e_k,
                mu=float(mu_guess),
                kBT=kbt0,
                delta=delta,
                spin_shift=1,
                dist_type=d0,
                U=float(meta["U"]),
            )
            plt.close("all")
            qc1, qc2 = st.columns(2)
            qw = float(qc1.slider("Quick view width", 4.0, 20.0, 7.0, 0.5))
            qh = float(qc2.slider("Quick view height", 3.0, 14.0, 5.0, 0.5))
            lim2d = st.checkbox("Quick view: custom color limits", value=False)
            vmin = vmax = None
            if lim2d:
                cl1, cl2 = st.columns(2)
                vmin = float(cl1.number_input("Color min", value=0.0))
                vmax = float(cl2.number_input("Color max", value=1.0))
            fig, ax = plt.subplots(figsize=(qw, qh))
            if meta["dim"] == "2D" and ky is not None:
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
            else:
                ax.plot(kx, np.diag(occ) if occ.ndim == 2 else occ, lw=2)
                ax.set_xlabel(r"$k$")
                ax.set_ylabel("occupancy")
            ax.set_title(str(d0))
            st.pyplot(fig, clear_figure=True)
        else:
            slice_type = st.selectbox("1D cut", ["horizontal", "diagonal"], index=0)
            n_opts = list(n_arr.astype(float))
            k_opts = list(kBT_arr.astype(float))
            n_pick = st.multiselect("n values", options=n_opts, default=n_opts[: min(3, len(n_opts))])
            k_pick = st.multiselect("k_BT values", options=k_opts, default=k_opts[: min(2, len(k_opts))])
            if not n_pick or not k_pick:
                st.error("Select at least one value for **n** and for **k_BT**.")
            else:
                H_plot = float(st.selectbox("H_a for plot", list(H_arr.astype(float))))
                n_kd_flag = st.checkbox("n_kd channel / correlator view (Hubbard tooling)", value=False)
                y0, y1 = st.slider("y axis (occupancy)", -0.05, 2.2, (-0.05, 1.15))
                fig_w = st.slider("Figure width (1D)", 10, 28, int(max(10, style["figsize"][0] * 1.8)))
                fig_h = st.slider("Figure height (1D)", 4, 16, int(max(4, style["figsize"][1])))
                fig2_w = st.slider("Figure width (2D)", 8, 28, int(max(8, style["figsize"][0] * 1.4)))
                fig2_h = st.slider("Figure height (2D)", 4, 18, int(max(4, style["figsize"][1])))

                names = [d.name if hasattr(d, "name") else str(d) for d in distros]
                plt.close("all")
                plot_distributions_notebook_style(
                    Plotter=PlotterCls,
                    e_k=e_k,
                    kx=kx,
                    ky=ky if ky is not None else kx,
                    dim=2 if meta["dim"] == "2D" else 1,
                    n_list=[float(x) for x in n_pick],
                    kBT_list=[float(x) for x in k_pick],
                    H=H_plot,
                    U=float(meta["U"]),
                    n_kd=n_kd_flag,
                    t_prim=float(meta["t_prim"]),
                    distributions=names,
                    slice_type=slice_type,
                    ylims=(float(y0), float(y1)),
                    figsize_1d=(float(fig_w), float(fig_h)),
                    figsize_2d=(float(fig2_w), float(fig2_h)),
                )
                show_current_figure(st)
                if not bool(style["text_on_plot"]):
                    st.caption(
                        f"Dystrybucje: n={len(n_pick)} punktow, kBT={len(k_pick)} punktow, "
                        f"H_a={H_plot:.6g}, slice={slice_type}."
                    )


if __name__ == "__main__":
    main()
