"""
Entropy calculation module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .distributions import distribution


def safe_log(x: np.ndarray) -> np.ndarray:
    """Compute the logarithm safely, avoiding log(0) issues."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(x > 0, np.log(x), 0.0)


def entropy_single_channel(n: np.ndarray, distro: str) -> np.ndarray:
    """
    Calculate entropy per k-state for a single spin channel.
    Parameters
    ----------
    n : np.ndarray
        Occupancy per spin channel (values in [0, 1]).
    distro : str
        Distribution type ("FD", "FLB", "NFLB", "NFLB_i").

    Returns
    -------
    np.ndarray
        Entropy values per k-state.
    """
    # Precompute logs
    ln_n = safe_log(n)
    ln_1_n = safe_log(1 - n)

    if distro == "FD":
        # Standard fermion entropy: -n ln n - (1-n) ln (1-n)
        return -n * ln_n - (1 - n) * ln_1_n  # type: ignore[no-any-return]

    if distro == "FLB":
        # FLB: S = S_fermion + n ln 2
        return -n * ln_n - (1 - n) * ln_1_n + n * np.log(2)  # type: ignore[no-any-return]

    if distro in ["NFLB", "NFLB_i"]:
        # NFLB (valid for n_tot <= 1, i.e., n_single <= 0.5)
        # S_single(x) = -x ln x + x ln 2 + (1-x) ln(1-x) - (1-2x) ln(1-2x)
        ln_1_2n = safe_log(1 - 2 * n)
        return -n * ln_n + n * np.log(2) + (1 - n) * ln_1_n - (1 - 2 * n) * ln_1_2n  # type: ignore[no-any-return]

    raise ValueError(f"Unknown distribution: {distro}")


def entropy(
    df: pd.DataFrame,
    E: np.ndarray,
    g: float = 2.0,
    muB_eV_per_T: float = 5.788381806e-5,
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

    entropies = []

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

        # Handle temperature
        beta = np.inf if k_BT == 0 else 1.0 / k_BT

        if H_a == 0:
            # Spin degenerate case
            # distribution() returns total occupancy (0..2 for FD)
            n_single = distribution(E, mu, beta, 0.0, 0, distro_str)

            # Total entropy = 2 * single_channel_entropy
            s_k = 2.0 * entropy_single_channel(n_single, distro_str)
            avg_entropy = np.mean(s_k)
        else:
            # Zeeman split case
            delta = 0.5 * g * muB_eV_per_T * H_a

            E_up = E - delta
            E_down = E + delta

            # distribution() returns total occupancy form (0..2)
            # Physical single-spin occupancy is distribution()/2
            n_up = distribution(E_up, mu, beta, 0.0, 1, distro_str)
            n_down = distribution(E_down, mu, beta, 0.0, -1, distro_str)

            # Total entropy = S(n_up) + S(n_down)
            s_k = entropy_single_channel(n_up, distro_str) + entropy_single_channel(
                n_down, distro_str
            )
            avg_entropy = np.mean(s_k)

        entropies.append(avg_entropy)

    # Return new DataFrame with entropy column
    df_out = df.copy()
    df_out["entropy"] = entropies
    return df_out
