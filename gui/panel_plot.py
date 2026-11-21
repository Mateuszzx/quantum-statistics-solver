from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from qss.distributions import DistType


def panel_plot(  # noqa: PLR0913, PLR0912, PLR0915
    df: pd.DataFrame,
    *,
    y_col: str,
    y_label: str,
    x_var: str,
    x_scale: str,
    panel_var: str | None,
    n_cols: int,
    fixed: dict[str, float],
    distros: list[DistType],
    style: dict[str, object],
) -> None:
    """
    Generic panel plot from the (distro, n, k_BT, H_a) multi-index results.

    - X axis: x_var
    - Panels: panel_var (if None -> single plot)
    - Lines: distributions (one line per selected distro)
    - Remaining variables are fixed via `fixed`
    """
    dfp = df.reset_index()
    if "distro" in dfp.columns:
        dfp["distro"] = dfp["distro"].astype(str).map(lambda s: s.split(".")[-1])
    for k, v in fixed.items():
        dfp = dfp[dfp[k] == float(v)]

    x_scale = str(x_scale).lower().strip()
    if x_scale not in ("linear", "log"):
        x_scale = "linear"
    if x_scale == "log" and (dfp[x_var].astype(float) <= 0).any():
        x_scale = "linear"

    latex_var: dict[str, str] = {
        "n": r"$n$",
        "mu": r"$\mu/|t|$",
        "k_BT": r"$k_B T/|t|$",
        "H_a": r"$H_a$ (T)",
        "U": r"$U/|t|$",
    }

    def _var_label(v: str) -> str:
        return latex_var.get(v, v)

    if panel_var is None:
        panel_vals: list[float | str] = ["__single__"]
    else:
        panel_vals = list(pd.unique(dfp[panel_var]))

    n_panels = len(panel_vals)
    n_cols = max(1, int(n_cols))
    n_rows = int(np.ceil(n_panels / n_cols))

    fig_w, fig_h = (float(style["figsize"][0]), float(style["figsize"][1]))
    xlims = style.get("xlims")
    ylims = style.get("ylims")

    plt.close("all")
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(fig_w * n_cols, fig_h * n_rows),
        squeeze=False,
        sharex=False,
        sharey=False,
    )

    distro_names = [d.name if hasattr(d, "name") else str(d).split(".")[-1] for d in distros]
    # Match `tools/scientific_plot_v2.py` defaults.
    style_map: dict[str, tuple[object, object]] = {
        "FD": ("blue", "-"),
        "FLB": ("red", "--"),
        "NFL": ("purple", "-."),
        "NFLB": ("green", ":"),
        "SSL": ("black", (0, (3, 1, 1, 1))),  # dotted dash
        "HUBBARD": ("orange", (0, (5, 1))),  # long dashes
        "HUBBARD ($n_{kd}$)": ("brown", (0, (1, 1))),
    }
    fallback_colors = [plt.cm.tab10(i % 10) for i in range(max(1, len(distro_names)))]
    fallback_map = {
        name: fallback_colors[i % len(fallback_colors)] for i, name in enumerate(distro_names)
    }

    def _pretty(v: object) -> str:
        if isinstance(v, (float, np.floating, int, np.integer)):
            return f"{float(v):.6g}"
        return str(v)

    legend_handles = None
    legend_labels = None
    for i, pv in enumerate(panel_vals):
        ax = axes[i // n_cols][i % n_cols]
        if panel_var is None:
            dfi = dfp
            title = ""
        else:
            dfi = dfp[dfp[panel_var] == pv]
            title = rf"{_var_label(str(panel_var))}={_pretty(pv)}"

        for dname in distro_names:
            dsub = dfi[dfi["distro"] == dname].sort_values(by=x_var)
            if dsub.empty:
                continue
            color, ls = style_map.get(dname, (fallback_map[dname], "-"))
            ax.plot(
                dsub[x_var].astype(float),
                dsub[y_col].astype(float),
                lw=2,
                label=dname,
                color=color,
                linestyle=ls,
            )

        ax.set_xlabel(_var_label(str(x_var)), fontsize=14)
        ax.set_ylabel(y_label, fontsize=14)
        if title:
            ax.set_title(title, fontsize=16)
        ax.grid(True, alpha=0.25)
        ax.set_xscale(x_scale)
        if xlims is not None:
            ax.set_xlim(tuple(xlims))  # type: ignore[arg-type]
        if ylims is not None:
            ax.set_ylim(tuple(ylims))  # type: ignore[arg-type]

        if legend_handles is None:
            legend_handles, legend_labels = ax.get_legend_handles_labels()

    for j in range(n_panels, n_rows * n_cols):
        axes[j // n_cols][j % n_cols].set_visible(False)

    if legend_handles and legend_labels:
        fig.legend(
            legend_handles,
            legend_labels,
            loc="lower center",
            ncol=min(len(legend_labels), max(1, n_cols * 2)),
            bbox_to_anchor=(0.5, -0.01),
            frameon=False,
            fontsize=18,
        )
        fig.tight_layout(rect=(0, 0.05, 1, 1))
    else:
        fig.tight_layout()
