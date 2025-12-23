"""
Distribution functions for quantum statistics.

This module provides convenient wrappers around the C++ _qss library
for computing quantum statistical distributions (Fermi-Dirac, etc.)
on numpy arrays.
"""

from __future__ import annotations

from enum import Enum
from typing import cast

import numpy as np
from numpy.typing import NDArray

EnergyData = NDArray[np.floating] | float | int
LatticeData = NDArray[np.floating] | float | int


try:
    from . import _qss  # type: ignore[attr-defined]
except ImportError as e:
    raise ImportError("The _qss extension module is not available.") from e


class DistType(Enum):
    """
    Quantum statistical distribution types.

    Attributes
    ----------
    FD : Fermi-Dirac distribution
    FLB : Fermi Liquid with Boltzmann distribution
    NFLB : Non-Fermi Liquid with Boltzmann distribution
    """

    FD = "FD"
    FLB = "FLB"
    NFL = "NFL"
    NFLB = "NFLB"
    SSL = "SSL"


def prepare_2D_lattice_data(
    kx: LatticeData,
    ky: LatticeData,
    t: float = -1.0,
    t_prim: float = 0.25,
) -> tuple[NDArray[np.floating], NDArray[np.floating], NDArray[np.floating]]:
    """
    Prepare 2D lattice data for a simple 2D tight-binding model.

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
    tuple[np.ndarray, np.ndarray, np.ndarray]
        `(e_k, KX, KY)` where `e_k` is the energy dispersion (same shape as the
        meshgrid), and `KX`, `KY` are the meshgrid arrays of wave vectors.
    """
    kx_array = np.atleast_1d(np.asarray(kx, dtype=float))
    ky_array = np.atleast_1d(np.asarray(ky, dtype=float))
    KX, KY = np.meshgrid(kx_array, ky_array)

    # 2D tight-binding dispersion relation
    cos_kx = np.cos(KX)
    cos_ky = np.cos(KY)
    e_k = 2 * t * (cos_kx + cos_ky) + 2 * t_prim * (cos_kx + cos_ky)

    return e_k, KX, KY


def _get_cpp_dist_type(dist_type: str | DistType | _qss.DistType) -> _qss.DistType:
    """Convert distribution type to C++ DistType."""

    if isinstance(dist_type, _qss.DistType):
        return dist_type

    if isinstance(dist_type, DistType):
        dist_type = dist_type.value

    # String conversion
    dist_map = {
        "FD": _qss.DistType.FD,
        "FLB": _qss.DistType.FLB,
        "NFLB": _qss.DistType.NFLB,
        "SSL": _qss.DistType.SSL,
        "NFL": _qss.DistType.NFL,
    }

    if dist_type not in dist_map:
        raise ValueError(
            f"Invalid distribution type: {dist_type}. " f"Must be one of {list(dist_map.keys())}"
        )

    return dist_map[dist_type]


def distribution(  # noqa: PLR0913
    e_k: np.ndarray,
    mu: float,
    beta: float,
    delta: float,
    spin_shift: int,
    dist_type: str | DistType,
) -> np.ndarray:
    """
    Compute distribution values for given energies using specified distribution type.

    Parameters
    ----------
    e_k : np.ndarray
        Energy values (can be 1D or 2D array)
    mu : float
        Chemical potential
    beta : float
        Inverse temperature (1/kT)
    delta : float
        Zeeman shift
    spin_shift : int
        Spin shift (+1 or -1)
    dist_type : str or DistType
        Distribution type: 'FD', 'FLB', 'NFLB' or DistType enum value

    Returns
    -------
    np.ndarray
        Distribution values with the same shape as e_k

    Raises
    ------
    RuntimeError
        If C++ extension is not available
    ValueError
        If invalid distribution type or array dimensions

    Examples
    --------
    >>> import numpy as np
    >>> from qss.distributions import distribution, DistType
    >>> energies = np.linspace(-2, 2, 100)
    >>> fd_dist = distribution(energies, mu=0.0, beta=2.0, delta=0.0, spin_shift=1, dist_type='FD')
    >>> # or using enum
    >>> fd_dist = distribution(energies, mu=0.0, beta=2.0, delta=0.0, spin_shift=1, dist_type=DistType.FD)
    >>>
    >>> # 2D example
    >>> energies_2d = np.random.randn(10, 20)
    >>> flb_dist = distribution(energies_2d, mu=0.5, beta=1.5, delta=0.1, spin_shift=1, dist_type=DistType.FLB)
    """
    e_k = np.asarray(e_k, dtype=np.float64)
    cpp_dist_type = _get_cpp_dist_type(dist_type)

    if e_k.ndim == 1:
        return cast(
            NDArray[np.floating], _qss.distribution(cpp_dist_type, mu, beta, delta, spin_shift, e_k)
        )
    if e_k.ndim == 2:
        return cast(
            NDArray[np.floating], _qss.distribution(cpp_dist_type, mu, beta, delta, spin_shift, e_k)
        )
    raise ValueError(f"e_k must be 1D or 2D array, got shape {e_k.shape}")


__all__ = [
    "DistType",
    "distribution",
    "prepare_2D_lattice_data",
]
