"""Data formatting utilities for quantum statistics calculations."""

from typing import Any, List, Optional

import numpy as np
from numpy.typing import NDArray
import pandas as pd


def numpy4d_to_dataframe(
    mu: NDArray[np.floating],
    distros: List[Any],
    n: NDArray[np.floating],
    beta: NDArray[np.floating],
    H_a: NDArray[np.floating],
    distros_labels: Optional[List[str]] = None,
) -> pd.DataFrame:
    """
    Convert 4D numpy array to pandas DataFrame with MultiIndex.

    Parameters
    ----------
    mu : ndarray
        4D array with shape (ND, NN, NB, NH)
    distros : list
        List of distribution types (length ND)
    n : ndarray
        Array of particle densities (length NN)
    beta : ndarray
        Array of inverse temperatures (length NB)
    H_a : ndarray
        Array of applied fields (length NH)
    distros_labels : list of str, optional
        Custom labels for distribution types

    Returns
    -------
    pd.DataFrame
        DataFrame with MultiIndex
    """
    if mu.ndim != 4:
        msg = f"mu must be a 4-D array, got {mu.ndim}-D"
        raise ValueError(msg)

    ND, NN, NB, NH = mu.shape
    if len(n) != NN or len(distros) != ND or len(beta) != NB or len(H_a) != NH:
        msg = (
            f"Lengths of distros, n, beta, H_a must match mu.shape. "
            f"Got ND={ND}, NN={NN}, NB={NB}, NH={NH}"
        )
        raise ValueError(msg)

    # Extract distribution names
    if distros_labels is not None:
        dist_names = distros_labels
    elif hasattr(distros[0], "name"):
        dist_names = [str(getattr(d, "name", d)) for d in distros]
    else:
        dist_names = [str(d) for d in distros]

    # Reshape mu
    mu_reshaped = mu.transpose(0, 2, 3, 1).reshape(-1, NN)

    # Create MultiIndex
    index = pd.MultiIndex.from_product(
        [dist_names, beta, H_a], names=["distributions", "beta", "H_a"]
    )

    return pd.DataFrame(mu_reshaped, index=index, columns=n)
