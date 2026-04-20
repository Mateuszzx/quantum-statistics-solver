"""
Entropy calculation module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .distributions import distribution


def n_kd(
    df: pd.DataFrame,
    E: np.ndarray,
    g: float = 2.0,
    muB_eV_per_T: float = 1.0,
    U: float = 0.0,
) -> pd.DataFrame:
    """
    Compute average entropy for each state in the chemical potential DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing chemical potential results.
        Expected index levels: 'distro', 'n', 'k_BT', 'H_a'.
        Expected columns: 'mu', ...
    E : np.ndarray
        Energy dispersion array.
    g : float
        Landé g-factor.
    muB_eV_per_T : float
        Bohr magneton in eV/T.

    Returns
    -------
    pd.DataFrame
        Original DataFrame with a new 'entropy' column.
    """
    # Ensure E is numpy array
    E = np.asarray(E)

    nkd = []

    # Iterate over the DataFrame
    for row in df.itertuples():
        # Extract parameters from index or columns
        try:
            if isinstance(row.Index, tuple) and len(row.Index) == 4:
                distro_str, _, k_BT, H_a = row.Index
            else:
                distro_str = row.distro
                k_BT = row.k_BT
                H_a = row.H_a
        except AttributeError as e:
            raise ValueError(
                "Input DataFrame must have 'distro', 'n', 'k_BT', 'H_a' either as MultiIndex levels or columns."
            ) from e

        mu = row.mu

        delta = 0.5 * g * muB_eV_per_T * H_a

        match distro_str:
            case "HUBBARD":
                n_kd = distribution(
                    E, mu, k_BT, delta=delta, spin_shift=-1, dist_type=distro_str, U=U, n_kd=True
                )
            case _:
                n_up = distribution(
                    E, mu, k_BT, delta=delta, spin_shift=1, dist_type=distro_str, U=U
                )
                n_down = distribution(
                    E, mu, k_BT, delta=delta, spin_shift=-1, dist_type=distro_str, U=U
                )
                n_kd = n_up * n_down

        avg_nkd = np.mean(n_kd)

        nkd.append(avg_nkd)

    # Return new DataFrame with entropy column
    df_out = df.copy()
    df_out["nkd"] = nkd
    return df_out
