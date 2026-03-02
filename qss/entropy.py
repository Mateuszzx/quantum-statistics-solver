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
        Distribution type ("FD", "FLB", "NFLB").

    Returns
    -------
    np.ndarray
        Entropy values per k-state.
    """
    # Precompute logs
    ln_n = safe_log(n)
    ln_1_n = safe_log(1 - n)

    match distro:
        case "FD":
            # Standard fermion entropy: -n ln n - (1-n) ln (1-n)
            return -n * ln_n - (1 - n) * ln_1_n  # type: ignore[no-any-return]
        case "FLB":
            # FLB: S = S_fermion + n ln 2
            return -n * ln_n - (1 - n) * ln_1_n + n * np.log(2)  # type: ignore[no-any-return]
        case "NFLB":
            # NFLB (valid for n_tot <= 1, i.e., n_single <= 0.5)
            # S_single(x) = -x ln x + x ln 2 + (1-x) ln(1-x) - (1-2x) ln(1-2x)
            ln_1_2n = safe_log(1 - 2 * n)
            return -n * ln_n + n * np.log(2) + (1 - n) * ln_1_n - (1 - 2 * n) * ln_1_2n  # type: ignore[no-any-return]
        case _:
            raise ValueError(f"Unknown distribution: {distro}")


def calculate_entropy_k(
    n_up: np.ndarray, n_down: np.ndarray, n_kd: np.ndarray, distro: str
) -> np.ndarray:
    """
    Calculate total entropy per k-state given occupancies of both channels.
    Handles correlations (like in SSL) correctly.

    Parameters
    ----------
    n_up : np.ndarray
        Occupancy for spin-up channel.
    n_down : np.ndarray
        Occupancy for spin-down channel.
    distro : str
        Distribution type.

    Returns
    -------
    np.ndarray
        Total entropy values per k-state.
    """
    match distro:
        case "NFL":
            # NFL: S = -n_up ln n_up - n_down ln n_down +
            #       2* (1 - n_up - n_down) ln(1 - n_up - n_down) +
            #       (1 - n_up) ln(1 - n_up)
            #       (1 - n_down) ln(1 - n_down)s

            n_hole = 1.0 - n_up - n_down
            # Ensure non-negative for log
            n_hole = np.maximum(n_hole, 0.0)
            ln_n_up = safe_log(n_up)
            ln_n_down = safe_log(n_down)
            ln_hole = safe_log(n_hole)
            ln_1_n_up = safe_log(1 - n_up)
            ln_1_n_down = safe_log(1 - n_down)

            return np.array(
                -n_up * ln_n_up
                - n_down * ln_n_down
                - 2 * n_hole * ln_hole
                + (1 - n_up) * ln_1_n_up
                + (1 - n_down) * ln_1_n_down
            )
        case "SSL":
            # SSL: S = -n_up ln n_up - n_down ln n_down - (1 - n_up - n_down) ln(1 - n_up - n_down)
            n_hole = 1.0 - n_up - n_down
            # Ensure non-negative for log
            n_hole = np.maximum(n_hole, 0.0)

            ln_n_up = safe_log(n_up)
            ln_n_down = safe_log(n_down)
            ln_hole = safe_log(n_hole)

            return np.array(-n_up * ln_n_up - n_down * ln_n_down - n_hole * ln_hole)

        case "HUBBARD":
            term1 = n_up - n_kd
            term2 = n_down - n_kd
            term3 = 1 - n_up - n_down + n_kd

            return np.array(
                -term1 * safe_log(term1)
                - term2 * safe_log(term2)
                - n_kd * safe_log(n_kd)
                - term3 * safe_log(term3)
            )
        case _:
            # Default: sum of independent channel entropies
            return np.array(
                entropy_single_channel(n_up, distro) + entropy_single_channel(n_down, distro)
            )


def entropy(
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

        if H_a == 0:
            # Spin degenerate case
            # distribution() returns single-spin occupancy (0..1 for FD)
            n_single = distribution(E, mu, k_BT, delta=0.0, spin_shift=1, dist_type=distro_str, U=U)
            n_kd = distribution(
                E, mu, k_BT, delta=0.0, spin_shift=-1, dist_type=distro_str, U=U, n_kd=True
            )
            # Total entropy
            # Pass same n for both channels (degenerate)
            s_k = calculate_entropy_k(n_single, n_single, n_kd, distro_str)
            avg_entropy = np.mean(s_k)
        else:
            # Zeeman split case
            delta = 0.5 * g * muB_eV_per_T * H_a

            # distribution() returns single-spin occupancy (0..1 for FD)
            n_up = distribution(E, mu, k_BT, delta=delta, spin_shift=1, dist_type=distro_str, U=U)
            n_down = distribution(
                E, mu, k_BT, delta=delta, spin_shift=-1, dist_type=distro_str, U=U
            )
            n_kd = distribution(
                E, mu, k_BT, delta=delta, spin_shift=-1, dist_type=distro_str, U=U, n_kd=True
            )
            # Total entropy
            s_k = calculate_entropy_k(n_up, n_down, n_kd, distro_str)
            avg_entropy = np.mean(s_k)

        entropies.append(avg_entropy)

    # Return new DataFrame with entropy column
    df_out = df.copy()
    df_out["entropy"] = entropies
    return df_out
