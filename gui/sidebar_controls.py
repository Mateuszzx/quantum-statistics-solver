"""Streamlit sidebar widgets shared across views."""

from __future__ import annotations

import streamlit as st

from gui.constants import DIM, DIST_MAP
from gui.defaults import cfg
from qss.distributions import DistType


def sidebar_distros(*, key_prefix: str = "") -> list[DistType]:
    st.sidebar.subheader("Distributions")
    names = list(DIST_MAP.keys())
    picked = st.sidebar.multiselect(
        "Select distributions",
        options=names,
        default=list(cfg("distributions", "default")),
        key=f"{key_prefix}distros",
    )
    return [DIST_MAP[k] for k in picked]


def sidebar_lattice(*, key_prefix: str = "") -> tuple[str, int, float, float]:
    st.sidebar.subheader("Lattice")
    dim = DIM
    n_cfg = cfg("lattice", "N")
    n_grid = int(
        st.sidebar.number_input(
            "N (k-grid per axis)",
            min_value=int(n_cfg["min"]),
            max_value=int(n_cfg["max"]),
            value=int(n_cfg["value"]),
            step=int(n_cfg["step"]),
            key=f"{key_prefix}N",
        )
    )
    t = float(st.sidebar.number_input("t (NN hopping)", value=float(cfg("lattice", "t")), key=f"{key_prefix}t"))
    t_prim = float(
        st.sidebar.number_input("t' (NNN hopping)", value=float(cfg("lattice", "t_prim")), key=f"{key_prefix}t_prim")
    )
    return dim, n_grid, t, t_prim
