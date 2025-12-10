""" Plot function for scientific articles/4letters"""
from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from qss import chemical_potential as cp
from qss import distributions as dist

DIST_COLORS = {
    "FD": "blue",
    "FLB": "red",
    "NFL": "purple",
    "NFLB": "green",
    "SSL": "black",
    "NFLB_i": "orange"
}


def _setup_axis(ax, text: Optional[str] = None, text_pos: Tuple[float, float] = (0.025, 0.985),
                xlabel: Optional[str] = None, ylabel: Optional[str] = None, fontsize: int = 12,
                ticks_len: int = 8, ticks_width: int = 3,
                legend: bool = True, legend_position: str = "best",
                xlims: Optional[Tuple[float, float]] = None,
                ylims: Optional[Tuple[float, float]] = None,
                show_log_vals: bool = False,
                xscale: str = "linear", yscale: str = "linear",
                spine_width: float = 1.5) -> None:
    """Helper function to setup axis styling."""

    if text:
        ax.text(text_pos[0], text_pos[1], text, ha='left', va='top', transform=ax.transAxes,
                fontdict={'fontsize': fontsize, 'color': 'black', 'fontweight': 'bold'})

    if xlabel:
        ax.set_xlabel(xlabel, fontsize=fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=fontsize)

    if legend:
        ax.legend(loc=legend_position, fontsize=fontsize)

    ax.grid(True)
    ax.tick_params(axis='both', which='major', direction="in", top=True, left=True,
                   right=True, bottom=True, length=ticks_len,
                   width=ticks_width, labelsize=fontsize)

    if xscale != "linear":
        ax.set_xscale(xscale)
    if yscale != "linear":
        ax.set_yscale(yscale)

    for spine in ax.spines.values():
        spine.set_linewidth(spine_width)

    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    if show_log_vals:
        ax2 = ax.twinx()
        log_vals = [np.log(2), 2*np.log(2), 3*np.log(2)]
        ax2.tick_params(axis='y', which='major', direction="in", width=ticks_width,
                        length=ticks_len, labelsize=fontsize)
        ax2.set_yticks(log_vals, labels=["$\\ln2$", "$2 \\ln2$", "$3 \\ln2$"], fontsize=fontsize)

        supp_colors = ['green', 'blue', 'red']
        # Assuming ylims are set on main ax, we sync them or use default for entropy
        # The original code hardcoded 3.2*ln(2)
        limit = 3.2 * np.log(2)
        ax.set_ylim(0, limit)
        ax2.set_ylim(0, limit)

        for i, log_v in enumerate(log_vals):
            ax2.axhline(log_v, color=supp_colors[i], linestyle='--', linewidth=1.5, alpha=0.7)


def plot_energy_distribution_along_k_const_n(e_k: np.ndarray, dim:int = 2, kx: Optional[np.ndarray] = None ,
                            ky: Optional[np.ndarray] = None , figsize: Tuple[int, int] = (15, 5),
                            n: np.ndarray = np.array([0.5, 0.9, 0.99]), N: int = 100,
                            betas: np.ndarray = np.array([1000, 100, 10]),
                            fontsize: int = 12, t_prim: float = 0.25,
                            ticks_len: int = 8, ticks_width: int = 3, ylims: Optional[Tuple[float, float]] = None,
                            text_pos: Tuple[float,float]=(0.025, .985),
                            legend: bool = True, legend_position: str = "best") -> None:


    n_num = len(n)
    num_betas = len(betas)
    _ , ax = plt.subplots(nrows=num_betas, ncols=n_num, figsize=(figsize[0] * n_num, figsize[1] * num_betas))
    distributions = ["FD", "FLB", "NFL", "NFLB", "SSL"]


    ax = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    if dim == 2 and kx is not None:
        for beta_i, beta in enumerate(betas):
           for j,n_i in enumerate(n):
                for distro in distributions:
                    # mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n_i, dim=dim, e_k=e_k, kx=kx, ky=ky)
                    df_mu = cp.chemical_potential(distros=[distro], n_vals=[n_i], beta_vals=[beta], H_a_vals=[0.0], E=e_k)
                    mu = df_mu["mu"].iloc[0]
                    n_k = dist.distribution(e_k, mu=mu, beta=beta, dist_type=distro)
                    linewidth = 2 if distro == "FD" else 1
                    color = DIST_COLORS.get(distro, "black")

                    ax[beta_i * n_num + j].plot(kx[N//2:], n_k[N//2,N//2:], color=color, marker="o", linewidth=linewidth, markersize=2, label=distro)

                # Set text
                text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$n = {n_i:.3f}$
$k_B T/|t| = {1/beta:.2f}$"""

                _setup_axis(ax[beta_i * n_num + j], text, text_pos, r"$\Gamma - X$", "occupation number",
                            fontsize, ticks_len, ticks_width, legend, legend_position, ylims=ylims)

                # Set proper tick labels for k_x
                ax[beta_i * n_num + j].set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])


        # plt.tight_layout()  # Adjust layout to make room for the suptitle
    else:
        raise ValueError("Invalid dimension or missing kx/ky for 2D plot.")



def plot_energy_distribution_diagonal_k_const_n(e_k: np.ndarray, dim:int = 2, kx: Optional[np.ndarray] = None ,
                            ky: Optional[np.ndarray] = None , figsize: Tuple[int, int] = (15, 5),
                            n: np.ndarray = np.array([0.5, 0.9, 0.99]), N: int = 100,
                            betas: np.ndarray = np.array([1000, 100, 10]),
                            fontsize: int = 12, t_prim: float = 0.25,
                            ticks_len: int = 8, ticks_width: int = 3, ylims: Optional[Tuple[float, float]] = None,
                            text_pos: Tuple[float,float]=(0.025, .985),
                            legend: bool = True, legend_position: str = "best") -> None:


    n_num = len(n)
    beta_num = len(betas)
    #####################################################################
    _ , ax = plt.subplots(nrows=beta_num, ncols=n_num, figsize=(figsize[0]*n_num, figsize[1]*beta_num))
    distributions = ["FD", "FLB", "NFL", "NFLB", "SSL"]

    ax = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    if dim == 2 and kx is not None:
        for beta_i, beta in enumerate(betas):
           for j,n_i in enumerate(n):
                for distro in distributions:
                    # mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n_i, dim=dim, e_k=e_k, kx=kx, ky=ky)
                    df_mu = cp.chemical_potential(distros=[distro], n_vals=[n_i], beta_vals=[beta], H_a_vals=[0.0], E=e_k)
                    mu = df_mu["mu"].iloc[0]
                    n_k = dist.distribution(e_k, mu=mu, beta=beta, dist_type=distro)
                    linewidth = 2 if distro == "FD" else 1
                    color = DIST_COLORS.get(distro, "black")

                    ax[beta_i*n_num + j].plot(kx[N//2:], np.diag(n_k)[N//2:], color=color, marker="o", linewidth=linewidth, markersize=2, label=distro)

                # Set text
                text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$n = {n_i:.3f}$
$k_B T/|t| = {1/beta:.2f}$"""

                _setup_axis(ax[beta_i * n_num + j], text, text_pos, r"$\Gamma - M$", "occupation number",
                            fontsize, ticks_len, ticks_width, legend, legend_position, ylims=ylims)

                # Set proper tick labels for k_x
                ax[beta_i*n_num + j].set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])


        # plt.tight_layout()  # Adjust layout to make room for the suptitle
    else:
        raise ValueError("Invalid dimension or missing kx/ky for 2D plot.")



def plot_energy_distribution(e_k: np.ndarray, dim:int = 2, kx: Optional[np.ndarray] = None ,
                            figsize: Tuple[int, int] = (10, 6), cmap: str = 'viridis',
                            ky: Optional[np.ndarray]= None, n: float = 0.8, beta: float = 1) -> None:

    vmax = [2,2,1]

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=figsize)
    distributions = ["FD", "FLB", "NFL", "NFLB", "SSL"]
    fig
    if dim == 2 and kx is not None and ky is not None:

        title = f"Statistical Distribution: $N = 100 \\times 100$, $n = {n:.2f}$,\n" +  "$ t = -1$, $t' / |t| = 0.25$," + " $\\frac{k_BT}{|t|} = $" + f"${1/beta:.2f}$"
        fig.suptitle(title, fontsize=16)
        for i, distro in enumerate(distributions):
            # mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n, dim=dim, e_k=e_k, kx=kx, ky=ky)
            df_mu = cp.chemical_potential(distros=[distro], n_vals=[n], beta_vals=[beta], H_a_vals=[0.0], E=e_k)
            mu = df_mu["mu"].iloc[0]
            n_k = dist.distribution(e_k, mu=mu, beta=beta, dist_type=distro)

            cf = ax[i].imshow(n_k, extent=[-np.pi, np.pi-0.01, -np.pi, np.pi-0.01], cmap=cmap, vmin=0, vmax=vmax[i])
            cf.set_clim(0, vmax[i])
            plt.colorbar(cf, ax=ax[i])

            ax[i].set_title(f"{distro}")
            ax[i].tick_params(axis='both', which='major', direction ="in", width = 0.01, labelsize=10)
            ax[i].set_xticks(np.array([-np.pi, -np.pi/2, 0, np.pi/2, np.pi-0.1]), ["$-\\pi$", "$-\\pi/2$", "0", "$\\pi/2$", "$\\pi$"])
            ax[i].set_yticks(np.array([-np.pi, -np.pi/2, 0, np.pi/2, np.pi-0.1]), ["$-\\pi$", "$-\\pi/2$", "0", "$\\pi/2$", "$\\pi$"])
            ax[i].set_xlabel("$k_x$")
            ax[i].set_ylabel("$k_y$")

        plt.tight_layout()  # Adjust layout to make room for the suptitle
    else:
        raise ValueError("Invalid dimension or missing kx/ky for 2D plot.")


def plot_result_n(results: pd.DataFrame, n: Union[List[float], np.ndarray],
                betas: Union[List[float], np.ndarray], y_label: str,
                distros: List[dist.DistType], H_a: float = 0.0, y_col: str = "mu", show_log_vals:bool = False,
                legend: bool = True, only_one_legend: bool = False, cols: int = 5,
                legend_position: str ="lower right", fontsize: int = 10,
                text_pos: Tuple[float,float]=(0.025, .985), xlims: Optional[Tuple[float, float]] = None,
                ylims: Optional[Tuple[float, float]] = None, t_prim: float = 0.25,
                ticks_len: int = 8, ticks_width: int=3, figsize: Tuple[int,int] = (12,6),
                ssl_theory:list = [], blue_line: bool = False) -> None:

    """Plot entropy as a function of density.

    Args:
        results (list): List of entropy values.
        n (list): List of density values.
        betas (list): List of beta values.
        distro (str): Distribution type.
        dim (int): Dimension (1 or 2).
        kx (Optional[np.ndarray]): Optional x-axis grid for 2D.
        ky (Optional[np.ndarray]): Optional y-axis grid for 2D.
"""
    # Prepare the figure
    rows = int(np.ceil(len(betas)/cols))

    _, ax = plt.subplots(nrows=rows, ncols=cols, figsize=(figsize[0]*cols, figsize[1]*rows))

    # Flatten axes for easy iteration
    ax_flat = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    for j,beta in enumerate(betas):
        k_BT = 1/beta
        # Plotting each distribution
        for distro in distros:
            d_label = distro.name if hasattr(distro, "name") else str(distro)
            try:
                # Select data for specific distro, k_BT, H_a
                subset = results.xs((d_label, k_BT, H_a), level=["distro", "k_BT", "H_a"])
                if "n" in subset.index.names:
                    subset = subset.sort_index(level="n")
                else:
                    subset = subset.sort_index()
                y_val = subset[y_col].values
            except KeyError:
                print(f"Warning: Data not found for {d_label}, beta={beta}, H_a={H_a}")
                continue

            y_val = y_val.reshape(-1)
            color = DIST_COLORS.get(d_label, "black")
            if blue_line and d_label == "FD":
                # Plotting the blue line for FD distribution
                ax_flat[j].plot(n, y_val, label=d_label, color=color, marker="o", markersize=4, linewidth=3)
            else:
                # Plotting the other distributions
                ax_flat[j].plot(n, y_val, label=d_label, color=color,  marker="o", markersize=3)

        if j in ssl_theory:
            # Plotting the SSL theory line
            n_vals = np.linspace(0.01, 0.99, 100)
            y_vals = n_vals * np.log(2) - n_vals * np.log(n_vals) - (1 - n_vals) * np.log(1 - n_vals)
            ax_flat[j].plot(n_vals, y_vals, color='black', linestyle='--', label="SSL", linewidth=1.5)

        # Set text
        text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$k_B T/|t| = {1/beta:.3f}$"""

        legend_arg = legend and (not only_one_legend or j == 0)
        _setup_axis(ax_flat[j], text, text_pos, "band filling, n", y_label,
                    fontsize, ticks_len, ticks_width, legend_arg, legend_position,
                    xlims, ylims, show_log_vals)

def plot_results_temperature(results: pd.DataFrame, n: Union[List[float], np.ndarray],
                betas: Union[List[float], np.ndarray], y_label: str,
                distros: List[dist.DistType], H_a: float = 0.0, y_col: str = "mu", show_log_vals:bool = False,
                legend: bool = True, only_one_legend: bool = False, cols: int = 5,
                legend_position = "lower right", figsize: Tuple[int,int]=(10,6),
                text_pos: Tuple[float,float]=(0.025, .985), xlims: Optional[Tuple[float, float]] = None,
                ylims: Optional[Tuple[float, float]] = None,t_prim: float = 0.25,
                fontsize:int = 10 ,ticks_len: int = 8, ticks_width: int=3, blue_line: bool = False) -> None:

    """Plot entropy as a function of density.

    Args:
        results (list): List of entropy values.
        n (list): List of density values.
        betas (list): List of beta values.
        distro (str): Distribution type.
        dim (int): Dimension (1 or 2).
        kx (Optional[np.ndarray]): Optional x-axis grid for 2D.
        ky (Optional[np.ndarray]): Optional y-axis grid for 2D.
    """
    # Prepare the figure
    plt.clf()
    rows = int(np.ceil(len(n)/cols))
    _ , ax = plt.subplots(nrows=rows, ncols=cols, figsize=(figsize[0]*cols, figsize[1]*rows))


    # Flatten axes for easy iteration
    ax_flat = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    for j, n_i in enumerate(n):

        # Plotting each distribution
        for distro in distros:
            d_label = distro.name if hasattr(distro, "name") else str(distro)

            try:
                subset = results.xs((d_label, n_i, H_a), level=["distro", "n", "H_a"])
                subset = subset.sort_index()
                y_plot = subset[y_col].values
                x_plot = subset.index.values
            except KeyError:
                print(f"Warning: Data not found for {d_label}, n={n_i}, H_a={H_a}")
                continue

            color = DIST_COLORS.get(d_label, "black")
            if blue_line and d_label == "FD":
                # Plotting the blue line for FD distribution
                ax_flat[j].semilogx(x_plot, y_plot , label=d_label, color=color, marker="o", markersize=4, linewidth=3)
            else:
                # Plotting the other distributions
                ax_flat[j].semilogx(x_plot, y_plot , label=d_label, color=color, marker="o", markersize=3)


        # Set text
        text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$ n = {n_i}$"""

        legend_arg = legend and (not only_one_legend or j == 0)
        _setup_axis(ax_flat[j], text, text_pos, "temperature $k_B  T / |t|$", y_label,
                    fontsize, ticks_len, ticks_width, legend_arg, legend_position,
                    xlims, ylims, show_log_vals, xscale="log")

def plot_results_Ha(results: pd.DataFrame, n: Union[List[float], np.ndarray],
                beta: float, y_label: str,
                distros: List[dist.DistType], y_col: str = "mu", show_log_vals:bool = False,
                legend: bool = True, only_one_legend: bool = False, cols: int = 5,
                legend_position = "lower right", figsize: Tuple[int,int]=(10,6),
                text_pos: Tuple[float,float]=(0.025, .985), xlims: Optional[Tuple[float, float]] = None,
                ylims: Optional[Tuple[float, float]] = None,t_prim: float = 0.25,
                fontsize:int = 10 ,ticks_len: int = 8, ticks_width: int=3, blue_line: bool = False) -> None:

    """Plot variable as a function of magnetic field H_a.

    Args:
        results (pd.DataFrame): Results DataFrame.
        n (list): List of density values (each creates a subplot).
        beta (float): Fixed beta value.
        y_label (str): Y-axis label.
        distros (list): List of distributions.
        y_col (str): Column to plot.
    """
    # Prepare the figure
    plt.clf()
    rows = int(np.ceil(len(n)/cols))
    _ , ax = plt.subplots(nrows=rows, ncols=cols, figsize=(figsize[0]*cols, figsize[1]*rows))

    # Flatten axes for easy iteration
    ax_flat = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    k_BT = 1.0 / beta

    for j, n_i in enumerate(n):

        # Plotting each distribution
        for distro in distros:
            d_label = distro.name if hasattr(distro, "name") else str(distro)

            try:
                # Select data for specific distro, n, k_BT
                # Index is (distro, n, k_BT, H_a)
                # We want to plot vs H_a
                subset = results.xs((d_label, n_i, k_BT), level=["distro", "n", "k_BT"])
                subset = subset.sort_index()
                y_plot = subset[y_col].values
                x_plot = subset.index.values # H_a values
            except KeyError:
                print(f"Warning: Data not found for {d_label}, n={n_i}, k_BT={k_BT}")
                continue

            color = DIST_COLORS.get(d_label, "black")
            if blue_line and d_label == "FD":
                # Plotting the blue line for FD distribution
                ax_flat[j].plot(x_plot, y_plot , label=d_label, color=color, marker="o", markersize=4, linewidth=3)
            else:
                # Plotting the other distributions
                ax_flat[j].plot(x_plot, y_plot , label=d_label, color=color, marker="o", markersize=3)


        # Set text
        text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$ n = {n_i}$
$k_B T/|t| = {k_BT:.3f}$"""

        legend_arg = legend and (not only_one_legend or j == 0)
        _setup_axis(ax_flat[j], text, text_pos, "magnetic field $H_a$", y_label,
                    fontsize, ticks_len, ticks_width, legend_arg, legend_position,
                    xlims, ylims, show_log_vals, xscale="log", spine_width=1.7)

def plot_3d_results(
    results: pd.DataFrame,
    row_param: str,
    row_vals: Union[List[float], np.ndarray],
    x_param: str,
    y_param: str,
    z_label: str,
    distros: List[Union[str, dist.DistType]] = ["FD", "FLB", "NFLB"],
    z_col: str = "mu",
    figsize: Tuple[int, int] = (15, 5),
    cmap: str = 'viridis',
    fontsize: int = 10
) -> None:
    """
    Create 3D surface plots of results.

    Args:
        results (pd.DataFrame): Results DataFrame with MultiIndex (distro, n, k_BT, H_a).
        row_param (str): Parameter to vary across rows ('n', 'k_BT', or 'H_a').
        row_vals (list): Values for the row parameter.
        x_param (str): Parameter for X axis ('n', 'k_BT', or 'H_a').
        y_param (str): Parameter for Y axis ('n', 'k_BT', or 'H_a').
        z_label (str): Label for Z axis.
        distros (list): List of distributions (columns).
        z_col (str): Column to plot on Z axis.
    """
    num_rows = len(row_vals)
    num_cols = len(distros)

    fig = plt.figure(figsize=(figsize[0], figsize[1] * num_rows))

    for i, r_val in enumerate(row_vals):
        for j, distro in enumerate(distros):
            d_label = distro.name if hasattr(distro, "name") else str(distro)

            try:
                # Filter by distro
                df_d = results.xs(d_label, level="distro", drop_level=False)

                # Filter by row_param value
                # Note: This requires exact float match.
                # If issues arise, we might need nearest match or tolerance.
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

                surf = ax.plot_surface(X, Y, Z, cmap=cmap, edgecolor='none', alpha=0.9)


                ax.set_xlabel(x_param, fontsize=fontsize)
                ax.set_ylabel(y_param, fontsize=fontsize)
                ax.set_zlabel(z_label, fontsize=fontsize)
                ax.set_title(f"{d_label}, {row_param}={r_val:.3g}", fontsize=fontsize)

            except Exception as e:
                print(f"Error plotting {d_label}, {row_param}={r_val}: {e}")
                # Add empty plot to maintain grid
                fig.add_subplot(num_rows, num_cols, i * num_cols + j + 1, projection='3d')

    plt.tight_layout()

def plot_3d_results_interactive(
    results: pd.DataFrame,
    row_param: str,
    row_vals: Union[List[float], np.ndarray],
    x_param: str,
    y_param: str,
    z_label: str,
    distros: List[Union[str, dist.DistType]] = ["FD", "FLB", "NFLB"],
    z_col: str = "mu",
    width: int = 1000,
    height: int = 400,
) -> None:
    """
    Create interactive 3D surface plots of results using Plotly.

    Args:
        results (pd.DataFrame): Results DataFrame with MultiIndex (distro, n, k_BT, H_a).
        row_param (str): Parameter to vary across rows ('n', 'k_BT', or 'H_a').
        row_vals (list): Values for the row parameter.
        x_param (str): Parameter for X axis ('n', 'k_BT', or 'H_a').
        y_param (str): Parameter for Y axis ('n', 'k_BT', or 'H_a').
        z_label (str): Label for Z axis.
        distros (list): List of distributions (columns).
        z_col (str): Column to plot on Z axis.
        width (int): Width of the plot.
        height (int): Height of each row.
    """
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
            d_label = distro.name if hasattr(distro, "name") else str(distro)
            subplot_titles.append(f"{d_label}, {row_param}={r_val:.3g}")

    fig = make_subplots(
        rows=num_rows, cols=num_cols,
        specs=[[{'type': 'surface'}] * num_cols] * num_rows,
        subplot_titles=subplot_titles,
        vertical_spacing=0.05,
        horizontal_spacing=0.02
    )

    for i, r_val in enumerate(row_vals):
        for j, distro in enumerate(distros):
            d_label = distro.name if hasattr(distro, "name") else str(distro)

            try:
                # Filter by distro
                df_d = results.xs(d_label, level="distro", drop_level=False)

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
                print(f"Error plotting {d_label}, {row_param}={r_val}: {e}")

    fig.update_layout(
        title_text="3D Results",
        height=height * num_rows,
        width=width,
    )
    fig.show()
