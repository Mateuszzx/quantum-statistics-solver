"""Plot layout controls for sweep panels."""

from __future__ import annotations

import streamlit as st

from gui.defaults import cfg


def plot_style_controls() -> dict[str, object]:
    p = cfg("plot")
    st.subheader("Plot customization")
    c1, c2, c3 = st.columns(3)
    fig_w = float(
        c1.number_input(
            "Figure width (per panel)",
            min_value=4.0,
            max_value=24.0,
            value=float(p["fig_width"]),
            step=0.5,
        )
    )
    fig_h = float(
        c2.number_input(
            "Figure height (per panel)",
            min_value=3.0,
            max_value=14.0,
            value=float(p["fig_height"]),
            step=0.5,
        )
    )
    cols = int(c3.number_input("Panels per row", min_value=1, max_value=8, value=int(p["cols"]), step=1))

    xl_enabled = st.checkbox("Set X axis limits", value=bool(p["x_limits_enabled"]))
    xlims: tuple[float, float] | None = None
    if xl_enabled:
        x1, x2 = st.columns(2)
        xmin = float(x1.number_input("x min", value=float(p["x_min"])))
        xmax = float(x2.number_input("x max", value=float(p["x_max"])))
        if xmax <= xmin:
            st.warning("x max should be greater than x min.")
        else:
            xlims = (xmin, xmax)

    yl_enabled = st.checkbox("Set Y axis limits", value=bool(p["y_limits_enabled"]))
    ylims: tuple[float, float] | None = None
    if yl_enabled:
        y1, y2 = st.columns(2)
        ymin = float(y1.number_input("y min", value=float(p["y_min"])))
        ymax = float(y2.number_input("y max", value=float(p["y_max"])))
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
