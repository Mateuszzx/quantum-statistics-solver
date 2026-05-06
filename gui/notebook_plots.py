"""Optional notebook-style figures (``ScientificPlotter`` from ``qss-results/tools``)."""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from gui.paths import ensure_qss_results_tools_on_path

if TYPE_CHECKING:
    import numpy as np
    from numpy.typing import NDArray


def load_scientific_plotter_v2():  # type: ignore[no-untyped-def]
    """Return ``ScientificPlotter`` class or ``None`` if the extra dependency layout is unavailable."""
    if ensure_qss_results_tools_on_path() is None:
        return None
    try:
        import tools.scientific_plot_v2 as sp  # type: ignore[import-untyped]

        return sp.ScientificPlotter
    except Exception:  # noqa: BLE001
        return None


def plot_distributions_notebook_style(  # noqa: PLR0913
    *,
    Plotter: Any,
    e_k: "NDArray[np.floating]",
    kx: "NDArray[np.floating]",
    ky: "NDArray[np.floating]",
    dim: int,
    n_list: list[float],
    kBT_list: list[float],
    H: float,
    U: float,
    n_kd: bool,
    t_prim: float,
    distributions: list[str],
    slice_type: str,
    ylims: tuple[float, float],
    figsize_1d: tuple[float, float],
    figsize_2d: tuple[float, float],
) -> None:
    """Call ``ScientificPlotter`` methods; requires matplotlib current figure context."""
    import matplotlib.pyplot as plt
    import numpy as np

    plotter = Plotter()
    if dim == 2 and len(distributions) == 1 and abs(H) > 1e-12:
        plotter.plot_energy_distribution_2d(
            e_k=e_k,
            dim=2,
            kx=kx,
            ky=ky,
            n=float(n_list[0]),
            kBT=float(kBT_list[0]),
            H=float(H),
            distributions=distributions,
            t_prim=float(t_prim),
            U=float(U),
            n_kd=bool(n_kd),
            N=int(e_k.shape[0]),
            figsize=figsize_2d,
        )
        return

    if dim == 2:
        plotter.plot_energy_distribution_1d(
            e_k=e_k,
            slice_type=slice_type,
            dim=2,
            kx=kx,
            n=np.asarray(n_list, dtype=float),
            kBTs=np.asarray(kBT_list, dtype=float),
            H=float(H),
            distributions=distributions,
            figsize=figsize_1d,
            ylims=ylims,
            U=float(U),
            n_kd=bool(n_kd),
            N=int(e_k.shape[0]),
        )
        return

    plt.figure()
    plt.text(0.1, 0.5, "Notebook-style 1D cuts require a 2D lattice.", fontsize=14)
    plt.axis("off")

