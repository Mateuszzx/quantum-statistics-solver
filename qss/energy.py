"""
Energy calculation module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .distributions import distribution
from .entropy import entropy as calculate_entropy


def energy(
    df: pd.DataFrame,
    E: np.ndarray,
    g: float = 2.0,
    muB_eV_per_T: float = 1.0,
    U: float = 0.0,
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
    U : float
        On-site interaction energy.

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
    for i, row in enumerate(df.itertuples()):
        # Extract parameters from index or columns
        try:
            if isinstance(row.Index, tuple):
                if len(row.Index) == 5:
                    distro_str, _, k_BT, H_a, row_U = row.Index
                elif len(row.Index) == 4:
                    distro_str, _, k_BT, H_a = row.Index
                    row_U = U
                else:
                    raise ValueError("Incorrect Index shape")
            else:
                distro_str = row.distro
                k_BT = row.k_BT
                H_a = row.H_a
                row_U = getattr(row, "U", U)
        except AttributeError as e:
            raise ValueError(
                "Input DataFrame must have 'distro', 'n', 'k_BT', 'H_a', 'U' either as MultiIndex levels or columns."
            ) from e

        mu = row.mu

        if H_a == 0:
            n_vals = distribution(E, mu, k_BT, delta=0.0, spin_shift=1, dist_type=distro_str, U=row_U)
            avg_energy = 2.0 * np.mean(E * n_vals)
        else:
            # Zeeman split case
            delta = 0.5 * g * muB_eV_per_T * H_a

            E_up = E - delta
            E_down = E + delta

            n_up = distribution(E, mu, k_BT, delta=delta, spin_shift=1, dist_type=distro_str, U=row_U)
            n_down = distribution(
                E, mu, k_BT, delta=delta, spin_shift=-1, dist_type=distro_str, U=row_U
            )
            avg_energy = np.mean(E_up * n_up + E_down * n_down)

        energies.append(avg_energy)

        # Calculate entropy if not present
        if hasattr(row, "entropy") and not pd.isna(row.entropy):
            entropy_val = row.entropy
        else:
            # We need to calculate entropy for this state.
            # Create a 1-row DataFrame preserving index
            # This is slow but handles the missing entropy case
            df_slice = df.iloc[[i]]
            # If entropy function adds 'entropy' column, use it.
            # Note: calculate_entropy requires mu which should be in df
            try:
                df_ent = calculate_entropy(df_slice, E, g, muB_eV_per_T, U)
                entropy_val = df_ent["entropy"].iloc[0]
            except Exception:
                # Fallback if entropy calculation fails (e.g. missing columns)
                entropy_val = 0.0

        free_energy.append(avg_energy - k_BT * entropy_val)

    # Return new DataFrame with energy column
    df_out = df.copy()
    df_out["energy"] = energies
    df_out["free_energy"] = free_energy
    return df_out
