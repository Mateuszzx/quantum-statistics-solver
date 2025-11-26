"""
Energy calculation module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .distributions import distribution
from .entropy import entropy as calculate_entropy, entropy_single_channel


def energy(
    df: pd.DataFrame,
    E: np.ndarray,
    g: float = 2.0,
    muB_eV_per_T: float = 5.788381806e-5,
) -> pd.DataFrame:
    """
    Compute average energy for each state in the chemical potential DataFrame.

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
        Original DataFrame with a new 'energy' column.
    """
    # Ensure E is numpy array
    E = np.asarray(E)
    
    if "entropy" not in df.columns:
        df = calculate_entropy(df, E, g=g, muB_eV_per_T=muB_eV_per_T)

    # We'll create a list of energies to append as a new column
    energies = []
    free_energy = []
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
        entropy = row.entropy

        # Handle temperature
        beta = np.inf if k_BT == 0 else 1.0 / k_BT

        if H_a == 0:
            # Spin degenerate case
            # distribution) returns total occupancy (factor of 2 included for FD)
            n_vals = distribution(E, mu, beta, distro_str)
            avg_energy = np.mean(E * n_vals)
        else:
            # Zeeman split case
            delta = 0.5 * g * muB_eV_per_T * H_a

            E_up = E - delta
            E_down = E + delta

            # distribution() returns 0..2 (for FD).
            # Physical occupancy per spin channel is distribution()/2.
            n_up = distribution(E_up, mu, beta, distro_str) / 2.0
            n_down = distribution(E_down, mu, beta, distro_str) / 2.0

            avg_energy = np.mean(E_up * n_up + E_down * n_down)

        energies.append(avg_energy)
        free_energy.append(avg_energy - k_BT * entropy)

    # Return new DataFrame with energy column
    df_out = df.copy()
    df_out["energy"] = energies
    df_out["free_energy"] = free_energy
    return df_out
