"""Distributions view: occupancy maps and 1D BZ cuts (no experiment run)."""

from __future__ import annotations

import matplotlib.pyplot as plt
import streamlit as st

from gui.constants import DIST_MAP
from gui.defaults import cfg
from gui.lattice import build_dispersion
from gui.mu_solver import solve_mu_from_n
from gui.occupancy_plot import plot_occupancy_map
from gui.scientific_plot import ScientificPlotter
from gui.sidebar_controls import sidebar_distros, sidebar_lattice
from gui.streamlit_utils import show_current_figure
from gui.sweep_grids import parse_float_list
from qss.chemical_potential import band_filling_from_mu


def view_distributions() -> None:  # noqa: PLR0912, PLR0915
    st.sidebar.caption("Standalone — no experiment run required.")
    dim, n_grid, t, t_prim = sidebar_lattice(key_prefix="dist_")
    distros = sidebar_distros(key_prefix="dist_")

    st.sidebar.subheader("State point")
    dv = cfg("distributions_view")
    prescribe_mu = st.sidebar.checkbox(
        r"Prescribe $\mu$ (skip $\mu(n)$ solve)",
        value=bool(dv["prescribe_mu"]),
        key="dist_prescribe_mu",
    )
    n_cfg = dv["n"]
    n_val: float | None = None
    mu_input: float | None = None
    if prescribe_mu:
        mu_input = float(
            st.sidebar.number_input(r"$\mu / |t|$", value=float(dv["mu"]), key="dist_mu")
        )
    else:
        n_val = float(
            st.sidebar.number_input(
                r"Band filling $n$",
                min_value=float(n_cfg["min"]),
                max_value=float(n_cfg["max"]),
                value=float(n_cfg["value"]),
                step=float(n_cfg["step"]),
                key="dist_n",
            )
        )
    kbt = float(st.sidebar.number_input(r"$k_B T / |t|$", value=float(dv["k_BT"]), key="dist_kbt"))
    h_a = float(st.sidebar.number_input(r"$H_a$ (T)", value=float(dv["H_a"]), key="dist_ha"))
    u_val = float(st.sidebar.number_input(r"$U / |t|$", value=float(dv["U"]), key="dist_u"))

    if not distros:
        st.info("Select at least one distribution in the sidebar.")
        return

    e_k, kx, ky = build_dispersion(dim, n_grid, t, t_prim)
    distro_names = [d.name if hasattr(d, "name") else str(d).split(".")[-1] for d in distros]

    st.header("Distributions")
    if prescribe_mu:
        st.caption(
            r"Occupancy from lattice + $(\mu, k_B T, H_a, U)$; "
            r"band filling $n$ is computed from occupancy. "
            "Does not use panel experiment results."
        )
    else:
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
    mu: float | None
    mu_err: str | None = None
    if prescribe_mu:
        assert mu_input is not None
        mu = float(mu_input)
        n_val = band_filling_from_mu(
            e_k, mu=mu, kBT=kbt, H_a=h_a, distro=str(d_pick), U=u_val
        )
    else:
        assert n_val is not None
        mu, mu_err = solve_mu_from_n(
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
        assert n_val is not None
        if prescribe_mu:
            st.caption(
                f"Prescribed μ / |t| = **{mu:.6g}** → n = {n_val:.6g}, "
                f"k_BT = {kbt:.6g}, H_a = {h_a:.6g}, U = {u_val:.6g}."
            )
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
        plot_occupancy_map(
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

    if dim == "2D" and ky is not None and not prescribe_mu:
        st.divider()
        st.subheader("1D BZ cuts")
        st.caption(
            "Uses `ScientificPlotter` — each curve fixes **n** and **$k_B T$** "
            "(μ solved per point). Independent of the Panels experiment."
        )
        slice_type = st.selectbox(
            "1D cut",
            ["horizontal", "diagonal"],
            index=int(dv["slice_index"]),
            key="dist_slice",
        )
        n_text = st.text_input(
            "n values (comma-separated)",
            value=str(dv["n_list"]),
            key="dist_n_list",
        )
        kbt_text = st.text_input(
            "$k_BT$ values (comma-separated)",
            value=str(dv["kbt_list"]),
            key="dist_kbt_list",
        )
        n_pick = parse_float_list(n_text) or parse_float_list(str(dv["n_list"]))
        k_pick = parse_float_list(kbt_text) or parse_float_list(str(dv["kbt_list"]))
        H_plot = float(st.number_input("$H_a$ for 1D plot", value=h_a, key="dist_h_plot"))
        U_plot = float(st.number_input("$U$ for 1D plot", value=u_val, key="dist_u_plot"))
        n_kd_flag = st.checkbox(
            "n_kd channel / correlator view (Hubbard tooling)",
            value=bool(dv["n_kd"]),
            key="dist_nkd",
        )
        ylim = tuple(float(x) for x in dv["ylim"])
        y0, y1 = st.slider("y axis (occupancy)", -0.05, 2.2, ylim, key="dist_ylim")
        fig_w = st.slider("Figure width (1D)", 10, 48, int(dv["fig_width"]), key="dist_fig_w")
        fig_h = st.slider("Figure height (1D)", 2, 8, int(dv["fig_height"]), key="dist_fig_h")

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
