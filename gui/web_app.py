"""
Streamlit GUI for quantum-statistics-solver: reproduces notebook workflows (sweeps in
``n``, ``k_B T``, ``H_a``; μ, entropy, energy, free energy, magnetization, ``n_{kd}``;
distribution / BZ cuts when ``qss-results/tools`` is available).

Run from the ``quantum-statistics-solver`` directory::

    streamlit run gui/web_app.py
"""

from __future__ import annotations

import streamlit as st

from gui.defaults import cfg
from gui.paths import quantum_statistics_solver_root
from gui.view_distributions import view_distributions
from gui.view_panels import view_panels

st.set_page_config(page_title="Quantum Statistics Solver", layout="wide")

SOLVER_ROOT = quantum_statistics_solver_root()


def main() -> None:
    st.title("Quantum statistics solver")
    st.caption(
        "Notebook-equivalent sweeps and BZ distributions — "
        f"project root: `{SOLVER_ROOT}`"
    )
    st.sidebar.caption(
        "Please cite: Kuterba, Wójcik & Spałek, arXiv:2611.26489 "
        "(CITATION.cff / README)."
    )
    st.sidebar.header("Quantum statistics solver")
    view = st.sidebar.radio(
        "View",
        ["Panels", "Distributions"],
        index=int(cfg("app", "default_view_index")),
        horizontal=True,
        key="main_view",
    )

    if view == "Panels":
        view_panels()
    else:
        view_distributions()


if __name__ == "__main__":
    main()
