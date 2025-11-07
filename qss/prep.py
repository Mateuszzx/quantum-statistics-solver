"""Distribution functions for quantum statistics."""

from typing import Tuple, Union

import numpy as np
from numpy.typing import NDArray

EnergyData = Union[NDArray[np.floating], float, int]
LatticeData = Union[NDArray[np.floating], float, int]


def prepare_2D_lattice_data(
    kx: LatticeData,
    ky: LatticeData,
    t: float = -1.0,
    t_prim: float = 0.25,
) -> Tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """
    Prepare 2D lattice data for tight-binding model.

    Parameters
    ----------
    kx : array-like or float
        x-component of wave vector
    ky : array-like or float
        y-component of wave vector
    t : float, optional
        Nearest-neighbor hopping parameter (default: -1.0)
    t_prim : float, optional
        Next-nearest-neighbor hopping parameter (default: 0.25)

    Returns
    -------
    tuple
        (e_k, KX, KY) where e_k is the energy dispersion,
        KX and KY are meshgrid arrays of wave vectors
    """
    kx_array: NDArray[np.floating] = (
        np.array(kx, dtype=float) if isinstance(kx, (int, float)) else np.asarray(kx, dtype=float)
    )
    ky_array: NDArray[np.floating] = (
        np.array(ky, dtype=float) if isinstance(ky, (int, float)) else np.asarray(ky, dtype=float)
    )
    KX, KY = np.meshgrid(kx_array, ky_array)

    # 2D tight-binding dispersion relation
    e_k: NDArray[np.floating] = 2 * t * (np.cos(KX) + np.cos(KY)) + 2 * t_prim * (
        np.cos(KX) + np.cos(KY)
    )

    return e_k, KX, KY
