""" Plot function for scientific articles/4letters"""
from typing import List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src import chemical_potential as cp
from src import distributions as dist


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
    distributions = dist.DISTRIBUTIONS[:3]
    col = ["blue", 'red', 'green']


    ax = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    if dim == 2 and kx is not None:
        for beta_i, beta in enumerate(betas):
           for j,n_i in enumerate(n):
                for i, distro in enumerate(distributions):
                    mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n_i, dim=dim, e_k=e_k, kx=kx, ky=ky)
                    n_k = dist.distribution(e_k, mu=mu, beta=beta, dist=distro)
                    linewidth = 2 if distro == "FD" else 1

                    ax[beta_i * n_num + j].plot(kx[N//2:], n_k[N//2,N//2:], color=col[i], marker="o", linewidth=linewidth, markersize=2, label=distro)

                # Set text
                text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$n = {n_i:.3f}$
$k_B T/|t| = {1/beta:.2f}$"""

                plt.text(text_pos[0], text_pos[1], text, ha='left', va='top', transform=ax[beta_i*n_num + j].transAxes,
                            fontdict={'fontsize': fontsize, 'color': 'black', 'fontweight': 'bold'})

                # Axis labels
                ax[beta_i * n_num + j].set_xlabel(r"$\Gamma - X$", fontsize=fontsize)
                ax[beta_i * n_num + j].set_ylabel("occupation number", fontsize=fontsize)
                # ax[beta_i * n_num + j].set_title(f"{distro}", fontsize=fontsize)

                # Legend
                if legend:
                    ax[beta_i * n_num + j].legend(loc=legend_position, fontsize=fontsize)

                # Grid & Ticks
                ax[beta_i * n_num + j].grid()
                ax[beta_i * n_num + j].tick_params(axis='both', which='major', direction ="in", top = True, left = True,
                                right = True, bottom = True, length =ticks_len,
                                width =ticks_width, labelsize=fontsize)

                # Set proper tick labels for k_x
                ax[beta_i * n_num + j].set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])

                # Spines
                for spine in ax[beta_i * n_num + j].spines.values():
                    spine.set_linewidth(1.5)

                if ylims:
                    ax[beta_i * n_num + j].set_ylim(ylims)


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
    distributions = dist.DISTRIBUTIONS[:3]
    col = ["blue", 'red', 'green']

    ax = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    if dim == 2 and kx is not None:
        for beta_i, beta in enumerate(betas):
           for j,n_i in enumerate(n):
                for i, distro in enumerate(distributions):
                    mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n_i, dim=dim, e_k=e_k, kx=kx, ky=ky)
                    n_k = dist.distribution(e_k, mu=mu, beta=beta, dist=distro)
                    linewidth = 2 if distro == "FD" else 1

                    ax[beta_i*n_num + j].plot(kx[N//2:], np.diag(n_k)[N//2:], color=col[i], marker="o", linewidth=linewidth, markersize=2, label=distro)

                # Set text
                text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$n = {n_i:.3f}$
$k_B T/|t| = {1/beta:.2f}$"""

                plt.text(text_pos[0], text_pos[1], text, ha='left', va='top', transform=ax[beta_i*n_num + j].transAxes,
                            fontdict={'fontsize': fontsize, 'color': 'black', 'fontweight': 'bold'})

                # Axis labels
                ax[beta_i*n_num + j].set_xlabel(r"$\Gamma - M$", fontsize=fontsize)
                ax[beta_i*n_num + j].set_ylabel("occupation number", fontsize=fontsize)
                # ax[beta_i*n_num + j].set_title(f"{distro}", fontsize=fontsize)

                # Legend
                if legend:
                    ax[beta_i*n_num + j].legend(loc=legend_position, fontsize=fontsize)

                # Grid & Ticks
                ax[beta_i*n_num + j].grid()
                ax[beta_i*n_num + j].tick_params(axis='both', which='major', direction ="in", top = True, left = True,
                                right = True, bottom = True, length =ticks_len,
                                width =ticks_width, labelsize=fontsize)

                # Set proper tick labels for k_x
                ax[beta_i*n_num + j].set_xticks(np.array([0, np.pi/2, np.pi]), ["0", "$\\pi/2$", "$\\pi$"])

                # Spines
                for spine in ax[beta_i*n_num + j].spines.values():
                    spine.set_linewidth(1.5)

                if ylims:
                    ax[beta_i*n_num + j].set_ylim(ylims)


        # plt.tight_layout()  # Adjust layout to make room for the suptitle
    else:
        raise ValueError("Invalid dimension or missing kx/ky for 2D plot.")



def plot_energy_distribution(e_k: np.ndarray, dim:int = 2, kx: Optional[np.ndarray] = None ,
                            figsize: Tuple[int, int] = (10, 6), cmap: str = 'viridis',
                            ky: Optional[np.ndarray]= None, n: float = 0.8, beta: float = 1) -> None:

    vmax = [2,2,1]

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=figsize)
    distributions = dist.DISTRIBUTIONS[:3]
    fig
    if dim == 2 and kx is not None and ky is not None:

        title = f"Statistical Distribution: $N = 100 \\times 100$, $n = {n:.2f}$,\n" +  "$ t = -1$, $t' / |t| = 0.25$," + " $\\frac{k_BT}{|t|} = $" + f"${1/beta:.2f}$"
        fig.suptitle(title, fontsize=16)
        for i, distro in enumerate(distributions):
            mu = cp.chemical_potential_computation(distro=distro, beta=beta, n=n, dim=dim, e_k=e_k, kx=kx, ky=ky)
            n_k = dist.distribution(e_k, mu=mu, beta=beta, dist=distro)

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
                distros: List[dist.DistType], show_log_vals:bool = False,
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
    col = ["blue", 'red', 'green']

    # Flatten axes for easy iteration
    ax_flat = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    for j,beta in enumerate(betas):
        # Plotting each distribution
        for i, distro in enumerate(distros):
            y_val = results.loc[(distro, beta), :]
            y_val = y_val.to_numpy()
            y_val = y_val.reshape(-1)
            if blue_line and distro == "FD":
                # Plotting the blue line for FD distribution
                ax_flat[j].plot(n, y_val, label=distro, color=col[i], marker="o", markersize=4, linewidth=3)
            else:
                # Plotting the other distributions
                ax_flat[j].plot(n, y_val, label=distro, color=col[i],  marker="o", markersize=3)

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

        plt.text(text_pos[0], text_pos[1], text, ha='left', va='top', transform=ax_flat[j].transAxes,
                    fontdict={'fontsize': fontsize, 'color': 'black', 'fontweight': 'bold'})

        # Axis labels
        ax_flat[j].set_xlabel("band filling, n", fontsize=fontsize)
        ax_flat[j].set_ylabel(y_label, fontsize=fontsize)

        # Legend
        if legend:
            if not only_one_legend:
                # ax_flat[j].legend()
                ax_flat[j].legend(loc=legend_position, fontsize=fontsize)
            else:
                ax_flat[0].legend(loc=legend_position, fontsize=fontsize)

        # Grid & Ticks
        ax_flat[j].grid()
        ax_flat[j].tick_params(axis='both', which='major', direction ="in", top = True, left = True,
                            right = True, bottom = True, length =ticks_len,
                            width =ticks_width,labelsize=fontsize)
        # Show log values
        if show_log_vals:
            ax2 = ax_flat[j].twinx()
            log_vals = [np.log(2), 2*np.log(2), 3*np.log(2)]
            ax2.tick_params(axis='y', which='major', direction ="in", width= ticks_width,
                            length = ticks_len, labelsize=fontsize)

            ax2.set_yticks(log_vals, labels=["$\\ln2$", "$2 ln2$", "$3 ln2$"], fontsize=fontsize)

            supp_colors = ['green', 'blue', 'red']
            ax_flat[j].set_ylim(0, 3.2*np.log(2))
            ax2.set_ylim(0, 3.2*np.log(2))
            for i,log in enumerate(log_vals):

                ax2.axhline(log, color=supp_colors[i], linestyle='--', linewidth=1.5, alpha=0.7)
        # Spines
        for spine in ax_flat[j].spines.values():
            spine.set_linewidth(1.5)

        # Xlims and Ylims
        if xlims is not None:
            ax_flat[j].set_xlim(xlims)

        if ylims is not None:
            ax_flat[j].set_ylim(ylims)

def plot_results_temperature(results: pd.DataFrame, n: Union[List[float], np.ndarray],
                betas: Union[List[float], np.ndarray], y_label: str,
                distros: List[dist.DistType], show_log_vals:bool = False,
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
    col = ["blue", 'red', 'green']

    # Flatten axes for easy iteration
    ax_flat = np.array(ax).flatten() if isinstance(ax, np.ndarray) else [ax]

    results.columns = results.columns.astype(float)
    # Invert Betas to get temperature
    x_val = np.array(betas[::-1])
    x_val = 1/x_val

    for j, n_i in enumerate(n):

        # Plotting each distribution
        for i, distro in enumerate(distros):
            y_val = results.loc[(distro, slice(None)),n_i]
            y_val = y_val[::-1]
            if blue_line and distro == "FD":
                # Plotting the blue line for FD distribution
                ax_flat[j].semilogx(x_val, y_val , label=distro, color=col[i], marker="o", markersize=4, linewidth=3)
            else:
                # Plotting the other distributions
                ax_flat[j].semilogx(x_val, y_val , label=distro, color=col[i], marker="o", markersize=3)


        # Set text
        text =f"""$t = -1$
$t' / |t| = {t_prim}$
$N = 100 \\times 100$
$ n = {n_i}$"""

        plt.text(text_pos[0], text_pos[1], text, ha='left', va='top', transform=ax_flat[j].transAxes,
                fontdict={'fontsize': fontsize, 'color': 'black', 'fontweight': 'bold'})

        # Axis labels
        ax_flat[j].set_ylabel(y_label, fontsize=fontsize)
        ax_flat[j].set_xlabel("temperature $k_B  T / |t|$", fontsize=fontsize)

        # Legend
        if legend:
            if not only_one_legend:
                ax_flat[j].legend(loc=legend_position, fontsize=fontsize)
            else:
                ax_flat[0].legend(loc=legend_position, fontsize=fontsize)

        # Grid & Ticks
        ax_flat[j].grid()
        ax_flat[j].tick_params(axis='both', which='major', direction ="in",  top = True, left = True,
                            right = True, bottom = True, width =ticks_width,
                            length = ticks_len, labelsize=fontsize)

        # Show log values
        if show_log_vals:
            ax2 = ax_flat[j].twinx()
            log_vals = [np.log(2), 2*np.log(2), 3*np.log(2)]
            ax2.tick_params(axis='y', which='major', direction ="in", width =ticks_width,
                            length = ticks_len, labelsize=fontsize)
            ax2.set_yticks(log_vals, labels=["$\\ln2$", "$2 ln2$", "$3 ln2$"], fontsize=fontsize)

            supp_colors = ['green', 'blue', 'red']
            ax_flat[j].set_ylim(0, 3.2*np.log(2))
            ax2.set_ylim(0, 3.2*np.log(2))
            for i,log in enumerate(log_vals):

                ax2.axhline(log, color=supp_colors[i], linestyle='--', linewidth=1.5, alpha=0.7)

        # Spines
        for spine in ax_flat[j].spines.values():
            spine.set_linewidth(1.5)

        # Xlims and Ylims
        if xlims is not None:
            ax_flat[j].set_xlim(xlims)

        if ylims is not None:
            ax_flat[j].set_ylim(ylims)
