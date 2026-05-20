"""Refactored plot function for scientific articles/letters"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np

from qss import chemical_potential
from qss import distributions as dist
from qss import solver_options as so


@dataclass
class PlotConfig:
    """Configuration for plot styling."""

    fontsize: int = 24
    ticks_len: int = 10
    ticks_width: int = 3
    spine_width: float = 3
    legend_fontsize: int = 16
    legend_loc: str = "best"
    cmap: str = "plasma"
    linewidth: float = 3.0

    dist_colors: Dict[str, str] = field(
        default_factory=lambda: {
            "FD": "blue",
            "FLB": "red",
            "NFL": "purple",
            "NFLB": "green",
            "SSL": "black",
            "HUBBARD": "orange",
            "HUBBARD ($n_{kd}$)": "brown",
        }
    )

    dist_linestyles: Dict[str, Union[str, Tuple]] = field(
        default_factory=lambda: {
            "FD": "-",
            "FLB": "--",
            "NFL": "-.",
            "NFLB": ":",
            "SSL": (0, (3, 1, 1, 1)),  # Dotted dash,
            "HUBBARD": (0, (5, 1)),  # Dashed with longer dashes
            "HUBBARD ($n_{kd}$)": (0, (1, 1)),
        }
    )

    def get_color(self, distro: str) -> str:
        return self.dist_colors.get(distro, "black")

    def get_linestyle(self, distro: str) -> Union[str, Tuple]:
        return self.dist_linestyles.get(distro, "-")


class ScientificPlotter:
    def __init__(self, config: Optional[PlotConfig] = None):
        self.config = config or PlotConfig()

    def _setup_axis(
        self,
        ax,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        title: Optional[str] = None,
        xlims: Optional[Tuple[float, float]] = None,
        ylims: Optional[Tuple[float, float]] = None,
        xscale: str = "linear",
        yscale: str = "linear",
        show_legend: bool = True,
        show_grid: bool = True,
    ) -> None:
        """Helper function to setup axis styling."""

        cfg = self.config

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

        ax.tick_params(
            axis="both",
            which="major",
            direction="in",
            top=True,
            left=True,
            right=True,
            bottom=True,
            length=cfg.ticks_len,
            width=cfg.ticks_width,
            labelsize=cfg.fontsize,
        )

        if xscale != "linear":
            ax.set_xscale(xscale)
        if yscale != "linear":
            ax.set_yscale(yscale)

        for spine in ax.spines.values():
            spine.set_linewidth(cfg.spine_width)
            spine.set_color("black")

        if xlims:
            ax.set_xlim(xlims)
        if ylims:
            ax.set_ylim(ylims)

    def _calculate_mu(
        self,
        e_k: np.ndarray,
        distro: str,
        n: float,
        kBT: float,
        H: float,
        g: float,
        mu_B: float,
        N: Optional[int] = None,
    ) -> float:
        _solver_options = (
            so.brent_defaults(
                abs_tol=1e-9,
                rel_tol=1e-6,
                expand_factor=1.4,
                initial_step=10.0,
                max_expands=64,
                max_iters=5000,
            )
        ).to_qss()

        df_mu = chemical_potential(
            distros=[distro],
            n_vals=[n],
            kBT_vals=[kBT],
            H_a_vals=[H],
            E=e_k,
            N=N,
            g=g,
            muB_eV_per_T=mu_B,
            solver_opts=_solver_options,
            x0=2.0,
            a0=-100.0,
            b0=100.0,
        )
        return float(df_mu["mu"].iloc[0])

    def add_bz_inset(
        self,
        ax,
        pos=(0.6, 0.55, 0.35, 0.35),
        highlight=[("G", "M")],
        highlight_color="red",
        highlight_lw=3.5,
        background_data=None,
    ) -> None:
        """Add a Brillouin-zone inset to an existing matplotlib Axes."""
        pts = {"G": (0, 0), "X": (1, 0), "M": (1, 1)}
        segs = [("G", "X"), ("X", "M"), ("M", "G"), ("G", "M")]

        # normalize highlight format
        highlight = [tuple(h.split("-")) if isinstance(h, str) else h for h in highlight]

        ax_inset = ax.inset_axes(pos)
        ax_inset.set_aspect("equal")
        ax_inset.set_xlim(-1, 1)
        ax_inset.set_ylim(-1, 1)
        ax_inset.set_xticks([])
        ax_inset.set_yticks([])

        if background_data is not None:
            ax_inset.imshow(
                background_data, extent=[-1, 1, -1, 1], origin="lower", cmap="plasma", alpha=0.8
            )

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
            "",
            xy=(1.3, 0.0),
            xytext=(-1.2, 0.0),
            arrowprops=dict(arrowstyle="-|>", linewidth=2.0, color="black"),
            zorder=10,
            clip_on=False,
            annotation_clip=False,
        )
        ax_inset.text(1.15, 0.05, r"$k_x$", ha="left", va="bottom", fontsize=14, clip_on=False)

        ax_inset.annotate(
            "",
            xy=(0.0, 1.3),
            xytext=(0.0, -1.2),
            arrowprops=dict(arrowstyle="-|>", linewidth=2.0, color="black"),
            zorder=10,
            clip_on=False,
            annotation_clip=False,
        )
        ax_inset.text(0.05, 1.1, r"$k_y$", ha="left", va="bottom", fontsize=14, clip_on=False)

        # --- Custom tick marks + labels on center axes ---
        tick_positions = [-1, 1]
        tick_labels = [r"$-\pi$", r"$\pi$"]
        tlen = 0.04
        for t, lbl in zip(tick_positions, tick_labels):
            ax_inset.plot([t, t], [-tlen, tlen], color="black", lw=1, zorder=9)
            ax_inset.text(t - 0.15, -0.05, lbl, ha="right", va="top", fontsize=14)
            ax_inset.plot([-tlen, tlen], [t, t], color="black", lw=1, zorder=9)
            ax_inset.text(-0.08, t - 0.1, lbl, ha="right", va="top", fontsize=14)

        return ax_inset

    def plot_energy_distribution_1d(
        self,
        e_k: np.ndarray,
        slice_type: str = "horizontal",  # "horizontal" or "diagonal"
        dim: int = 2,
        kx: Optional[np.ndarray] = None,
        n: Union[float, np.ndarray] = 0.9,
        kBTs: Union[float, np.ndarray] = 0.001,
        H: float = 0.0,
        g: float = 2.0,
        mu_B: float = 1,
        distributions: List[str] = ["FD", "FLB", "NFL"],
        figsize: Tuple[int, int] = (15, 5),
        ylims: Optional[Tuple[float, float]] = (-0.1, 1.1),
        bz_inset_pos: Tuple[float, float, float, float] = (0.65, 0.68, 0.25, 0.25),
        U: float = 0.0,
        n_kd: bool = False,
        N: Optional[int] = None,
    ) -> None:
        if dim != 2 or kx is None:
            raise ValueError("Invalid dimension or missing kx for 2D plot.")

        if N is None:
            N = e_k.shape[0]

        # Normalize inputs to arrays
        n_arr = np.array([n]) if np.isscalar(n) else np.array(n)
        kBTs_arr = np.array([kBTs]) if np.isscalar(kBTs) else np.array(kBTs)

        n_num = len(n_arr)
        num_kBTs = len(kBTs_arr)
        is_H_zero = np.isclose(H, 0.0)

        total_rows = n_num * (1 if is_H_zero else 2)

        _, ax = plt.subplots(
            nrows=num_kBTs, ncols=total_rows, figsize=(figsize[0], figsize[1] * total_rows)
        )
        ax = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

        delta = H * 0.5 * g * mu_B

        # Determine slice parameters
        if slice_type == "horizontal":
            slice_label = r"$\Gamma - X$"
            highlight = ["G-X"]
            # Assuming kx is symmetric and centered at 0, we want positive half
            # Original code: kx[N//2:]
            x_axis = kx[N // 2 :]
            slice_func = lambda arr: arr[N // 2, N // 2 :]
        elif slice_type == "diagonal":
            slice_label = r"$\Gamma - M$"
            highlight = ["G-M"]
            x_axis = kx[N // 2 :]  # Diagonal projection on kx is same length if square
            slice_func = lambda arr: np.diag(arr)[N // 2 :]
        else:
            raise ValueError(f"Unknown slice_type: {slice_type}")

        plot_idx = 0
        for kBT in kBTs_arr:
            for n_val in n_arr:
                # Prepare data for this configuration
                data_per_distro = {}
                for distro in distributions:
                    mu = self._calculate_mu(e_k, distro, n_val, kBT, H, g, mu_B, N)

                    n_k_plus = dist.distribution(
                        e_k,
                        mu=mu,
                        kBT=kBT,
                        delta=delta,
                        spin_shift=-1,
                        dist_type=distro,
                        U=U,
                        n_kd=n_kd,
                    )
                    n_k_minus = dist.distribution(
                        e_k,
                        mu=mu,
                        kBT=kBT,
                        delta=delta,
                        spin_shift=1,
                        dist_type=distro,
                        U=U,
                        n_kd=n_kd,
                    )

                    if distro in ["NFL", "NFLB", "SSL"]:
                        n_k_plus = np.clip(n_k_plus, 0, 1.0)
                        n_k_minus = np.clip(n_k_minus, 0, 1.0)

                    data_per_distro[distro] = (n_k_plus, n_k_minus)

                # Plotting

                if is_H_zero:
                    current_ax = ax[plot_idx]
                    for distro, (nkp, nkm) in data_per_distro.items():
                        total = nkp + nkm
                        y_vals = slice_func(total)
                        current_ax.plot(
                            x_axis,
                            y_vals,
                            color=self.config.get_color(distro),
                            linestyle=self.config.get_linestyle(distro),
                            linewidth=self.config.linewidth,
                            label=distro,
                        )

                    self.add_bz_inset(current_ax, pos=bz_inset_pos, highlight=highlight)
                    self._setup_axis(
                        current_ax, xlabel=slice_label, ylabel="spin occupancy $n_k$", ylims=ylims
                    )
                    current_ax.set_xticks(
                        np.array([0, np.pi / 2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"]
                    )
                    current_ax.set_xlim(0, np.pi)
                    plot_idx += 1
                else:
                    # Plot n_k+
                    ax_plus = ax[plot_idx]
                    for distro, (nkp, _) in data_per_distro.items():
                        y_vals = slice_func(nkp)
                        ax_plus.plot(
                            x_axis,
                            y_vals,
                            color=self.config.get_color(distro),
                            linestyle=self.config.get_linestyle(distro),
                            linewidth=self.config.linewidth,
                            label=distro,
                        )

                    self.add_bz_inset(ax_plus, pos=bz_inset_pos, highlight=highlight)
                    self._setup_axis(
                        ax_plus,
                        xlabel=slice_label,
                        ylabel="spin occupancy $n_{k\\uparrow}$",
                        ylims=ylims,
                    )
                    ax_plus.set_xticks(np.array([0, np.pi / 2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])
                    ax_plus.set_xlim(0, np.pi)
                    plot_idx += 1

                    # Plot n_k-
                    ax_minus = ax[plot_idx]
                    for distro, (_, nkm) in data_per_distro.items():
                        y_vals = slice_func(nkm)
                        ax_minus.plot(
                            x_axis,
                            y_vals,
                            color=self.config.get_color(distro),
                            linestyle=self.config.get_linestyle(distro),
                            linewidth=self.config.linewidth,
                            label=distro,
                        )

                    self.add_bz_inset(ax_minus, pos=bz_inset_pos, highlight=highlight)
                    self._setup_axis(
                        ax_minus,
                        xlabel=slice_label,
                        ylabel="spin occupancy $n_{k\\downarrow}$",
                        ylims=ylims,
                    )
                    ax_minus.set_xticks(
                        np.array([0, np.pi / 2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"]
                    )
                    ax_minus.set_xlim(0, np.pi)
                    plot_idx += 1
