"""Streamlit sidebar widgets shared across views."""

from __future__ import annotations

import streamlit as st

from gui.constants import DIM, DIST_MAP
from qss.distributions import DistType


def sidebar_distros(*, key_prefix: str = "") -> list[DistType]:
    st.sidebar.subheader("Distributions")
    names = list(DIST_MAP.keys())
    picked = st.sidebar.multiselect(
        "Select distributions",
        options=names,
        default=["FD", "NFL", "SSL"],
        key=f"{key_prefix}distros",
    )
    return [DIST_MAP[k] for k in picked]


def sidebar_lattice(*, key_prefix: str = "") -> tuple[str, int, float, float]:
    st.sidebar.subheader("Lattice")
    dim = DIM
    n_grid = int(
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
    return dim, n_grid, t, t_prim
