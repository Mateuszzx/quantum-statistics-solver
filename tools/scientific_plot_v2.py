"""Refactored plot function for scientific articles/letters"""
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union, Dict

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

from qss import chemical_potential
from qss import solver_options as so
from qss import distributions as dist

@dataclass
class PlotConfig:
    """Configuration for plot styling."""
    fontsize: int = 12
    ticks_len: int = 8
    ticks_width: int = 3
    spine_width: float = 1.5
    legend_fontsize: int = 12
    legend_loc: str = "best"
    cmap: str = "plasma"

    dist_colors: Dict[str, str] = field(default_factory=lambda: {
        "FD": "blue",
        "FLB": "red",
        "NFL": "purple",
        "NFLB": "green",
        "SSL": "black"
    })

    dist_linestyles: Dict[str, Union[str, Tuple]] = field(default_factory=lambda: {
        "FD": "-",
        "FLB": "--",
        "NFL": "-.",
        "NFLB": ":",
        "SSL": (0, (3, 1, 1, 1)) # Dotted dash
    })

    def get_color(self, distro: str) -> str:
        return self.dist_colors.get(distro, "black")

    def get_linestyle(self, distro: str) -> Union[str, Tuple]:
        return self.dist_linestyles.get(distro, "-")

class ScientificPlotter:
    def __init__(self, config: Optional[PlotConfig] = None):
        self.config = config or PlotConfig()

    def _setup_axis(self, ax,
                    xlabel: Optional[str] = None,
                    ylabel: Optional[str] = None,
                    title: Optional[str] = None,
                    text: Optional[str] = None,
                    text_pos: Tuple[float, float] = (0.025, 0.985),
                    xlims: Optional[Tuple[float, float]] = None,
                    ylims: Optional[Tuple[float, float]] = None,
                    xscale: str = "linear",
                    yscale: str = "linear",
                    show_legend: bool = True,
                    show_grid: bool = True,
                    show_log_vals: bool = False) -> None:
        """Helper function to setup axis styling."""

        cfg = self.config

        if text:
            ax.text(text_pos[0], text_pos[1], text, ha='left', va='top', transform=ax.transAxes,
                    fontdict={'fontsize': cfg.fontsize, 'color': 'black', 'fontweight': 'bold'})

        if xlabel:
            ax.set_xlabel(xlabel, fontsize=cfg.fontsize)
        if ylabel:
            ax.set_ylabel(ylabel, fontsize=cfg.fontsize)
        if title:
            ax.set_title(title, fontsize=cfg.fontsize)

        if show_legend:
            ax.legend(loc=cfg.legend_loc, fontsize=cfg.legend_fontsize)

        if show_grid:
            ax.grid(True)

        ax.tick_params(axis='both', which='major', direction="in", top=True, left=True,
                        right=True, bottom=True, length=cfg.ticks_len,
                        width=cfg.ticks_width, labelsize=cfg.fontsize)

        if xscale != "linear":
            ax.set_xscale(xscale)
        if yscale != "linear":
            ax.set_yscale(yscale)

        for spine in ax.spines.values():
            spine.set_linewidth(cfg.spine_width)

        if xlims:
            ax.set_xlim(xlims)
        if ylims:
            ax.set_ylim(ylims)

        if show_log_vals:
            self._add_log_vals_axis(ax)

    def _add_log_vals_axis(self, ax):
        cfg = self.config
        ax2 = ax.twinx()
        log_vals = [np.log(2), 2*np.log(2), 3*np.log(2)]
        ax2.tick_params(axis='y', which='major', direction="in", width=cfg.ticks_width,
                        length=cfg.ticks_len, labelsize=cfg.fontsize)
        ax2.set_yticks(log_vals, labels=["$\\ln2$", "$2 \\ln2$", "$3 \\ln2$"], fontsize=cfg.fontsize)

        supp_colors = ['green', 'blue', 'red']
        limit = 3.2 * np.log(2)
        # Sync ylims if possible, otherwise set default
        current_ylim = ax.get_ylim()
        if current_ylim[1] > limit or current_ylim[1] == 1.0: # Heuristic check
             ax.set_ylim(0, limit)

        ax2.set_ylim(ax.get_ylim())

        for i, log_v in enumerate(log_vals):
            ax2.axhline(log_v, color=supp_colors[i], linestyle='--', linewidth=1.5, alpha=0.7)

    def _calculate_mu(self, e_k: np.ndarray, distro: str, n: float, beta: float, H: float, g: float, mu_B: float, N: Optional[int] = None) -> float:
        _solver_options = (so.brent_defaults(abs_tol=1e-9,
                        rel_tol=1e-6,
                        expand_factor=1.4,
                        initial_step=10.0,
                        max_expands=64,
                        max_iters=5000 )).to_qss()

        df_mu = chemical_potential(
            distros=[distro],
            n_vals=[n],
            beta_vals=[beta],
            H_a_vals=[H],
            E=e_k,
            N=N,
            g=g,
            muB_eV_per_T=mu_B,
            solver_opts=_solver_options,
            x0=2.0,
            a0=-100.0,
            b0=100.0
        )
        return float(df_mu["mu"].iloc[0])

    def add_bz_inset(self, ax, pos=(0.6, 0.55, 0.35, 0.35), highlight=[("G","M")],
                     highlight_color="red", highlight_lw=3.5, background_data=None) -> None:
        """Add a Brillouin-zone inset to an existing matplotlib Axes."""
        pts = {"G": (0, 0), "X": (1, 0), "M": (1, 1)}
        segs = [("G", "X"), ("X", "M"), ("M", "G"), ("G", "M")]

        # normalize highlight format
        highlight = [tuple(h.split("-")) if isinstance(h, str) else h for h in highlight]

        ax_inset = ax.inset_axes(pos)
        ax_inset.set_aspect("equal")
        ax_inset.set_xlim(-1, 1)
        ax_inset.set_ylim(-1, 1)
        ax_inset.set_xticks([]); ax_inset.set_yticks([])

        if background_data is not None:
            ax_inset.imshow(background_data, extent=[-1, 1, -1, 1], origin='lower', cmap='plasma', alpha=0.8)

        # --- Draw segments ---
        for a, b in segs:
            x = [pts[a][0], pts[b][0]]
            y = [pts[a][1], pts[b][1]]
            if (a, b) in highlight or (b, a) in highlight:
                ax_inset.plot(x, y, color=highlight_color, lw=highlight_lw, zorder=2)
            else:
                ax_inset.plot(x, y, color="black", lw=1, zorder=1)

        # --- Labels for Γ, X, M ---
        label_pos = {"G": (-0.2, 0.2), "X": (0.8, 0.2), "M": (0.5, 0.8)}
        for name, (x, y) in pts.items():
            label = r"$\Gamma$" if name == "G" else f"${name}$"
            dx, dy = label_pos[name]
            ax_inset.text(dx, dy, label, fontsize=14, ha="center", va="center")

        # --- Center axes arrows (crossing borders) ---
        ax_inset.annotate(
            "", xy=(1.3, 0.0), xytext=(-1.2, 0.0),
            arrowprops=dict(arrowstyle="-|>", linewidth=2.0, color="black"),
            zorder=10, clip_on=False, annotation_clip=False
        )
        ax_inset.text(1.15, 0.05, r"$k_x$", ha="left", va="bottom", fontsize=14, clip_on=False)

        ax_inset.annotate(
            "", xy=(0.0, 1.3), xytext=(0.0, -1.2),
            arrowprops=dict(arrowstyle="-|>", linewidth=2.0, color="black"),
            zorder=10, clip_on=False, annotation_clip=False
        )
        ax_inset.text(0.05, 1.1, r"$k_y$", ha="left", va="bottom", fontsize=14, clip_on=False)

        # --- Custom tick marks + labels on center axes ---
        tick_positions = [-1, 1]
        tick_labels = [r"$-\pi$", r"$\pi$"]
        tlen = 0.04
        for t, lbl in zip(tick_positions, tick_labels):
            ax_inset.plot([t, t], [-tlen, tlen], color="black", lw=1, zorder=9)
            ax_inset.text(t-0.15, -0.05, lbl, ha="right", va="top", fontsize=14)
            ax_inset.plot([-tlen, tlen], [t, t], color="black", lw=1, zorder=9)
            ax_inset.text(-0.08, t-0.1, lbl, ha="right", va="top", fontsize=14)

        return ax_inset

    def plot_energy_distribution_1d(self, e_k: np.ndarray,
                                    slice_type: str = "horizontal", # "horizontal" or "diagonal"
                                    dim: int = 2,
                                    kx: Optional[np.ndarray] = None,
                                    n: Union[float, np.ndarray] = 0.9,
                                    betas: Union[float, np.ndarray] = 1000,
                                    H: float = 0.0, g: float = 2.0, mu_B: float = 1,
                                    distributions: List[str] = ["FD", "FLB", "NFL"],
                                    figsize: Tuple[int, int] = (15, 5),
                                    ylims: Optional[Tuple[float, float]] = (-0.1, 1.1),
                                    text_pos: Tuple[float,float]=(0.025, .985),
                                    bz_inset_pos: Tuple[float, float, float, float] = (0.65, 0.68, 0.25, 0.25),
                                    N: Optional[int] = None) -> None:

        if dim != 2 or kx is None:
             raise ValueError("Invalid dimension or missing kx for 2D plot.")

        if N is None:
            N = e_k.shape[0]

        # Normalize inputs to arrays
        n_arr = np.array([n]) if np.isscalar(n) else np.array(n)
        betas_arr = np.array([betas]) if np.isscalar(betas) else np.array(betas)

        n_num = len(n_arr)
        num_betas = len(betas_arr)
        is_H_zero = np.isclose(H, 0.0)

        total_rows = num_betas * n_num * (1 if is_H_zero else 2)

        fig, ax = plt.subplots(nrows=total_rows, ncols=1, figsize=(figsize[0], figsize[1] * total_rows))
        ax = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

        delta = H * 0.5 * g * mu_B

        # Determine slice parameters
        if slice_type == "horizontal":
            slice_label = r"$\Gamma - X$"
            highlight = ["G-X"]
            # Assuming kx is symmetric and centered at 0, we want positive half
            # Original code: kx[N//2:]
            x_axis = kx[N//2:]
            slice_func = lambda arr: arr[N//2, N//2:]
        elif slice_type == "diagonal":
            slice_label = r"$\Gamma - M$"
            highlight = ["G-M"]
            x_axis = kx[N//2:] # Diagonal projection on kx is same length if square
            slice_func = lambda arr: np.diag(arr)[N//2:]
        else:
            raise ValueError(f"Unknown slice_type: {slice_type}")

        plot_idx = 0
        for beta in betas_arr:
            for n_val in n_arr:
                # Prepare data for this configuration
                data_per_distro = {}
                for distro in distributions:
                    mu = self._calculate_mu(e_k, distro, n_val, beta, H, g, mu_B, N)

                    n_k_plus = dist.distribution(e_k, mu=mu, beta=beta, delta=delta, spin_shift=-1, dist_type=distro)
                    n_k_minus = dist.distribution(e_k, mu=mu, beta=beta, delta=delta, spin_shift=1, dist_type=distro)

                    if distro in ["NFL", "NFLB", "SSL"]:
                        n_k_plus = np.clip(n_k_plus, 0, 1.0)
                        n_k_minus = np.clip(n_k_minus, 0, 1.0)

                    data_per_distro[distro] = (n_k_plus, n_k_minus)

                # Plotting
                text = f"$n = {n_val:.2f}$\n$k_B T/|t| = {1/beta:.1e}$\n$\\mu H / |t| = {H}$"

                if is_H_zero:
                    current_ax = ax[plot_idx]
                    for distro, (nkp, nkm) in data_per_distro.items():
                        total = nkp + nkm
                        y_vals = slice_func(total)
                        current_ax.plot(x_axis, y_vals,
                                        color=self.config.get_color(distro),
                                        linestyle=self.config.get_linestyle(distro),
                                        linewidth=2, label=distro)

                    self.add_bz_inset(current_ax, pos=bz_inset_pos, highlight=highlight)
                    self._setup_axis(current_ax, xlabel=slice_label, ylabel="spin occupancy $n_k$",
                                     text=text, text_pos=text_pos, ylims=ylims)
                    current_ax.set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])
                    current_ax.set_xlim(0, np.pi)
                    plot_idx += 1
                else:
                    # Plot n_k+
                    ax_plus = ax[plot_idx]
                    for distro, (nkp, _) in data_per_distro.items():
                        y_vals = slice_func(nkp)
                        ax_plus.plot(x_axis, y_vals,
                                     color=self.config.get_color(distro),
                                     linestyle=self.config.get_linestyle(distro),
                                     linewidth=2, label=distro)

                    self.add_bz_inset(ax_plus, pos=bz_inset_pos, highlight=highlight)
                    self._setup_axis(ax_plus, xlabel=slice_label, ylabel="spin occupancy $n_{k+}$",
                                     text=text, text_pos=text_pos, ylims=ylims)
                    ax_plus.set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])
                    ax_plus.set_xlim(0, np.pi)
                    plot_idx += 1

                    # Plot n_k-
                    ax_minus = ax[plot_idx]
                    for distro, (_, nkm) in data_per_distro.items():
                        y_vals = slice_func(nkm)
                        ax_minus.plot(x_axis, y_vals,
                                      color=self.config.get_color(distro),
                                      linestyle=self.config.get_linestyle(distro),
                                      linewidth=2, label=distro)

                    self.add_bz_inset(ax_minus, pos=bz_inset_pos, highlight=highlight)
                    self._setup_axis(ax_minus, xlabel=slice_label, ylabel="spin occupancy $n_{k-}$",
                                     text=text, text_pos=text_pos, ylims=ylims)
                    ax_minus.set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])
                    ax_minus.set_xlim(0, np.pi)
                    plot_idx += 1

    def plot_energy_distribution_2d(self, e_k: np.ndarray,
                                    dim: int = 2,
                                    kx: Optional[np.ndarray] = None,
                                    ky: Optional[np.ndarray] = None,
                                    n: float = 0.8,
                                    beta: float = 1000,
                                    H: float = 0.0, g: float = 2.0, mu_B: float = 1,
                                    distributions: List[str] = ["FD", "FLB", "NFL"],
                                    figsize: Tuple[int, int] = (14, 6),
                                    t_prim: float = 0.25,
                                    N: Optional[int] = None) -> None:

        if dim != 2 or kx is None or ky is None:
            raise ValueError("Invalid dimension or missing kx/ky for 2D plot.")

        if N is None:
            N = e_k.shape[0]

        delta = H * 0.5 * g * mu_B
        extent = [kx.min(), np.pi, ky.min(), np.pi]
        ticks = [-np.pi, -np.pi/2, 0, np.pi/2, np.pi]
        tick_labels = ["$-\\pi$", "$-\\pi/2$", "0", "$\\pi/2$", "$\\pi$"]

        # If H is non-zero, we might want to plot spin up/down for a single distribution
        # or multiple distributions. The original code had two separate functions.
        # Let's handle the case where we plot multiple distributions side-by-side (H=0 usually)
        # OR single distribution spin-split (H!=0).

        if not np.isclose(H, 0.0) and len(distributions) == 1:
            # Spin split plot for single distribution
            distro = distributions[0]
            mu = self._calculate_mu(e_k, distro, n, beta, H, g, mu_B, N)
            n_k_plus = dist.distribution(e_k, mu=mu, beta=beta, delta=delta, spin_shift=-1, dist_type=distro)
            n_k_minus = dist.distribution(e_k, mu=mu, beta=beta, delta=delta, spin_shift=1, dist_type=distro)

            fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=figsize, sharey=True)

            im1 = ax1.imshow(n_k_plus, aspect='auto', extent=extent, origin='lower', cmap=self.config.cmap, vmin=0, vmax=1)
            self._setup_axis(ax1, xlabel="$k_x$", ylabel="$k_y$", title=f"$n_{{k+}}$, $H = {H}$ T",
                             show_legend=False, show_grid=False)

            im2 = ax2.imshow(n_k_minus, aspect='auto', extent=extent, origin='lower', cmap=self.config.cmap, vmin=0, vmax=1)
            self._setup_axis(ax2, xlabel="$k_x$", title=f"$n_{{k-}}$, $H = {H}$ T",
                             show_legend=False, show_grid=False)

            for ax in [ax1, ax2]:
                ax.set_xticks(ticks)
                ax.set_xticklabels(tick_labels)
                ax.set_yticks(ticks)
                ax.set_yticklabels(tick_labels)

            cbar = fig.colorbar(im2, ax=[ax1, ax2], label='spin occupancy', fraction=0.046, pad=0.04)
            cbar.ax.tick_params(labelsize=self.config.fontsize)

            title = f"Distribution: {distro}, $n = {n:.3f}$, $k_B T/|t| = {1/beta:.2f}$, $t'/|t| = {t_prim}$"
            fig.suptitle(title, fontsize=self.config.fontsize+2)

        else:
            # Multiple distributions (usually H=0, but can be used for H!=0 summing spins)
            fig, ax = plt.subplots(nrows=1, ncols=len(distributions), figsize=figsize)
            if len(distributions) == 1:
                ax = [ax]

            vmax_vals = [2, 2, 1] # Heuristic from original code

            for i, distro in enumerate(distributions):
                mu = self._calculate_mu(e_k, distro, n, beta, H, g, mu_B, N)
                # Summing spins if H=0 or just total occupancy
                n_k = dist.distribution(e_k, mu=mu, beta=beta, delta=delta, spin_shift=-1, dist_type=distro) + \
                      dist.distribution(e_k, mu=mu, beta=beta, delta=delta, spin_shift=1, dist_type=distro)

                vm = vmax_vals[i] if i < len(vmax_vals) else 1.0
                cf = ax[i].imshow(n_k, extent=extent, origin='lower', aspect='auto', cmap=self.config.cmap, vmin=0, vmax=vm)

                self._setup_axis(ax[i], xlabel="$k_x$", ylabel="$k_y$" if i==0 else None, title=f"{distro}",
                                 show_legend=False, show_grid=False)

                ax[i].set_xticks(ticks)
                ax[i].set_xticklabels(tick_labels)
                ax[i].set_yticks(ticks)
                ax[i].set_yticklabels(tick_labels)

                plt.colorbar(cf, ax=ax[i])

            title = f"Statistical Distribution: $N = {N} \\times {N}$, $n = {n:.2f}$,\n" + \
                    f"$ t = -1$, $t' / |t| = {t_prim}$," + f" $\\frac{{k_BT}}{{|t|}} = {1/beta:.2f}$"
            fig.suptitle(title, fontsize=16)
            plt.tight_layout()

    def plot_scalar_result(self, results: pd.DataFrame,
                           x_param: str, # "n", "k_BT", "H_a"
                           y_col: str,
                           y_label: str,
                           fixed_params: Dict[str, Union[float, int]], # e.g. {"beta": 100, "H_a": 0}
                           distros: List[str],
                           subplots_param: Optional[str] = None, # Parameter to create subplots for (e.g. "n" or "beta")
                           subplots_vals: Optional[List[float]] = None,
                           cols: int = 5,
                           figsize: Tuple[int, int] = (12, 6),
                           xlims: Optional[Tuple[float, float]] = None,
                           ylims: Optional[Tuple[float, float]] = None,
                           text_pos: Tuple[float, float] = (0.025, 0.985),
                           show_log_vals: bool = False,
                           xscale: str = "linear",
                           blue_line: bool = False,
                           ssl_theory: bool = False,
                           t_prim: float = 0.25,
                           N: int = 100) -> None:

        if subplots_param and subplots_vals is not None:
            rows = int(np.ceil(len(subplots_vals)/cols))
            fig, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(figsize[0]*cols, figsize[1]*rows))
            ax_flat = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]
        else:
            fig, ax = plt.subplots(figsize=figsize)
            ax_flat = [ax]
            subplots_vals = [None] # Dummy iteration

        for j, sp_val in enumerate(subplots_vals):
            current_ax = ax_flat[j]

            # Build context for filtering
            current_params = fixed_params.copy()
            if subplots_param:
                current_params[subplots_param] = sp_val

            # Extract needed params for text label
            n_val = current_params.get("n", "var")
            beta_val = current_params.get("beta", "var")
            H_val = current_params.get("H_a", 0.0)
            k_BT_val = 1.0/beta_val if isinstance(beta_val, (int, float)) else "var"

            for distro in distros:
                try:
                    # Filter logic depends on index structure. Assuming MultiIndex (distro, n, k_BT, H_a)
                    # We need to construct the slicer

                    # This part is tricky because the original code had specific logic for each plot type.
                    # We'll try to be generic.

                    # Create a slicer dictionary
                    slicer = {"distro": distro}
                    if "n" in current_params: slicer["n"] = current_params["n"]
                    if "beta" in current_params: slicer["k_BT"] = 1.0/current_params["beta"]
                    if "H_a" in current_params: slicer["H_a"] = current_params["H_a"]

                    # We need to handle floating point matching for index
                    # This is a simplified approach, might need robust float matching
                    subset = results.xs(distro, level="distro")

                    if "n" in slicer: subset = subset.xs(slicer["n"], level="n")
                    if "k_BT" in slicer: subset = subset.xs(slicer["k_BT"], level="k_BT")
                    if "H_a" in slicer: subset = subset.xs(slicer["H_a"], level="H_a")

                    subset = subset.sort_index()
                    y_vals = subset[y_col].values
                    x_vals = subset.index.values

                    color = self.config.get_color(distro)
                    if blue_line and distro == "FD":
                         current_ax.plot(x_vals, y_vals, label=distro, color=color, linewidth=2)
                    else:
                         current_ax.plot(x_vals, y_vals, label=distro, color=color, linewidth=2)

                except KeyError:
                    # print(f"Warning: Data not found for {distro} with params {current_params}")
                    continue

            if ssl_theory and x_param == "n":
                 n_vals = np.linspace(0.01, 0.99, 100)
                 y_theory = n_vals * np.log(2) - n_vals * np.log(n_vals) - (1 - n_vals) * np.log(1 - n_vals)
                 current_ax.plot(n_vals, y_theory, color='black', linestyle='--', label="SSL", linewidth=1.5)

            # Text label
            text = f"$t = -1$\n$t' / |t| = {t_prim}$\n$N = {N} \\times {N}$"
            if isinstance(n_val, (int, float)): text += f"\n$n = {n_val}$"
            if isinstance(k_BT_val, (int, float)): text += f"\n$k_B T/|t| = {k_BT_val:.2e}$"
            if isinstance(H_val, (int, float)) and H_val != 0: text += f"\n$H_a = {H_val}$"

            xlabel_map = {
                "n": "band filling, n",
                "k_BT": "temperature $k_B T / |t|$",
                "H_a": "magnetic field $H_a$"
            }

            self._setup_axis(current_ax, xlabel=xlabel_map.get(x_param, x_param), ylabel=y_label,
                             text=text, text_pos=text_pos, xlims=xlims, ylims=ylims,
                             show_log_vals=show_log_vals, xscale=xscale,
                             show_legend=(j==0))

    def plot_3d_results(self,
        results: pd.DataFrame,
        row_param: str,
        row_vals: Union[List[float], np.ndarray],
        x_param: str,
        y_param: str,
        z_label: str,
        distros: List[str] = ["FD", "FLB", "NFLB"],
        z_col: str = "mu",
        figsize: Tuple[int, int] = (15, 5),
        fontsize: int = 10
    ) -> None:
        """Create 3D surface plots of results."""
        num_rows = len(row_vals)
        num_cols = len(distros)

        fig = plt.figure(figsize=(figsize[0], figsize[1] * num_rows))

        for i, r_val in enumerate(row_vals):
            for j, distro in enumerate(distros):
                try:
                    # Filter by distro
                    df_d = results.xs(distro, level="distro", drop_level=False)
                    # Filter by row_param value
                    df_r = df_d.xs(r_val, level=row_param, drop_level=False)
                    # Reset index to make pivoting easy
                    data = df_r.reset_index()
                    # Pivot
                    pivot = data.pivot(index=y_param, columns=x_param, values=z_col)

                    X_vals = pivot.columns.values
                    Y_vals = pivot.index.values
                    X, Y = np.meshgrid(X_vals, Y_vals)
                    Z = pivot.values

                    # Add subplot
                    ax = fig.add_subplot(num_rows, num_cols, i * num_cols + j + 1, projection='3d')
                    ax.plot_surface(X, Y, Z, cmap=self.config.cmap, edgecolor='none', alpha=0.9)

                    ax.set_xlabel(x_param, fontsize=fontsize)
                    ax.set_ylabel(y_param, fontsize=fontsize)
                    ax.set_zlabel(z_label, fontsize=fontsize)
                    ax.set_title(f"{distro}, {row_param}={r_val:.3g}", fontsize=fontsize)

                except Exception as e:
                    print(f"Error plotting {distro}, {row_param}={r_val}: {e}")
                    # Add empty plot to maintain grid
                    fig.add_subplot(num_rows, num_cols, i * num_cols + j + 1, projection='3d')

        plt.tight_layout()

    def plot_3d_results_interactive(self,
        results: pd.DataFrame,
        row_param: str,
        row_vals: Union[List[float], np.ndarray],
        x_param: str,
        y_param: str,
        z_label: str,
        distros: List[str] = ["FD", "FLB", "NFLB"],
        z_col: str = "mu",
        width: int = 1000,
        height: int = 400,
    ) -> None:
        """Create interactive 3D surface plots of results using Plotly."""
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            print("Plotly is required for interactive plots. Please install it with 'pip install plotly'.")
            return

        num_rows = len(row_vals)
        num_cols = len(distros)

        # Create subplot titles
        subplot_titles = []
        for r_val in row_vals:
            for distro in distros:
                subplot_titles.append(f"{distro}, {row_param}={r_val:.3g}")

        fig = make_subplots(
            rows=num_rows, cols=num_cols,
            specs=[[{'type': 'surface'}] * num_cols] * num_rows,
            subplot_titles=subplot_titles,
            vertical_spacing=0.05,
            horizontal_spacing=0.02
        )

        for i, r_val in enumerate(row_vals):
            for j, distro in enumerate(distros):
                try:
                    # Filter by distro
                    df_d = results.xs(distro, level="distro", drop_level=False)
                    # Filter by row_param value
                    df_r = df_d.xs(r_val, level=row_param, drop_level=False)
                    # Reset index to make pivoting easy
                    data = df_r.reset_index()
                    # Pivot
                    pivot = data.pivot(index=y_param, columns=x_param, values=z_col)

                    X = pivot.columns.values
                    Y = pivot.index.values
                    Z = pivot.values

                    # Add surface trace
                    fig.add_trace(
                        go.Surface(x=X, y=Y, z=Z, colorscale='Viridis', showscale=False),
                        row=i + 1, col=j + 1
                    )

                    # Update scene labels
                    scene_num = i * num_cols + j + 1
                    scene_name = 'scene' if scene_num == 1 else f'scene{scene_num}'

                    fig.update_layout({
                        scene_name: {
                            'xaxis_title': x_param,
                            'yaxis_title': y_param,
                            'zaxis_title': z_label,
                        }
                    })

                except Exception as e:
                    print(f"Error plotting {distro}, {row_param}={r_val}: {e}")

        fig.update_layout(
            title_text="3D Results",
            height=height * num_rows,
            width=width,
        )
        fig.show()
