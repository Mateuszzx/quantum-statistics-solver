"""Occupancy map and 1D cut plots for the Distributions view."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import streamlit as st

from gui.constants import DIM_2D
from qss.distributions import DistType, distribution


def plot_occupancy_map(  # noqa: PLR0913
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
