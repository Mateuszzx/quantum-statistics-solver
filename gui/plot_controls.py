"""Plot layout controls for sweep panels."""

from __future__ import annotations

import streamlit as st


def plot_style_controls() -> dict[str, object]:
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
