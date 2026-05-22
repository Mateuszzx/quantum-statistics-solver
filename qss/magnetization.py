import numpy as np
import pandas as pd

from .distributions import distribution


def magnetization_single(n_up: np.ndarray, n_down: np.ndarray) -> float:
    if n_up.shape != n_down.shape:
        raise ValueError("n_up and n_down must have the same shape.")

    kx, ky = n_up.shape
    return float(np.abs(np.sum(n_up - n_down) / (kx * ky)))


def calculate_magnetization(
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

    Returns
    -------
    pd.DataFrame
        Original DataFrame with a new 'magnetization' column.
    """
    # Ensure E is numpy array
    E = np.asarray(E)

    # We'll create a list of energies to append as a new column
    magnetization = []

    # Iterate over the DataFrame
    for row in df.itertuples():
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
            mag = 0.0
        else:
            # Zeeman split case
            delta = 0.5 * g * muB_eV_per_T * H_a

            n_up = distribution(E, mu, k_BT, delta=delta, spin_shift=1, dist_type=distro_str, U=row_U)
            n_down = distribution(
                E, mu, k_BT, delta=delta, spin_shift=-1, dist_type=distro_str, U=row_U
            )
            # Incorrect formulation for Energy?
            mag = magnetization_single(n_up, n_down)

        magnetization.append(mag)
    # Return new DataFrame with energy column
    df_out = df.copy()
    df_out["magnetization"] = magnetization
    return df_out
