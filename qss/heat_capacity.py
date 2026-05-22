"""
Heat capacity calculation module.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .entropy import entropy as calculate_entropy


def heat_capacity(
    df: pd.DataFrame,
    E: np.ndarray,
    g: float = 2.0,
    muB_eV_per_T: float = 1.0,
    U: float = 0.0,
) -> pd.DataFrame:
    """
    Compute specific heat capacity for each state in the DataFrame.
    C = T * (\partial S / \partial T) = k_B T * (\partial S / \partial k_B T)

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing chemical potential and entropy results.
        Expected index levels: 'distro', 'n', 'k_BT', 'H_a'.
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
        Original DataFrame with a new 'heat_capacity' column.
    """
    if "entropy" not in df.columns:
        df = calculate_entropy(df, E, g=g, muB_eV_per_T=muB_eV_per_T)

    # We need to sort by the parameters to properly compute the derivative
    # Usually we take derivative w.r.t k_BT keeping n, H_a, distro, U constant
    
    # Store original row order/index
    original_idx = df.index
    
    # Reset index to easily group and sort (if index is MultiIndex)
    is_multiindex = isinstance(df.index, pd.MultiIndex)
    df_reset = df.reset_index()

    if "U" not in df_reset.columns:
        df_reset["U"] = U

    # Group by the ones that should be constant
    group_cols = ["distro", "n", "H_a"]
    if "U" in df_reset.columns:
        group_cols.append("U")
        
    # Calculate heat capacity
    # C = k_BT * dS / d(k_BT)
    heat_cap_list = np.zeros(len(df_reset))
    
    for name, group in df_reset.groupby(group_cols):
        # Sort by k_BT to compute derivative
        sorted_group = group.sort_values("k_BT")
        indices = sorted_group.index
        k_BT_vals = sorted_group["k_BT"].values
        S_vals = sorted_group["entropy"].values
        
        if len(k_BT_vals) > 1:
            # np.gradient handles non-uniform spacing
            dS_dkBT = np.gradient(S_vals, k_BT_vals)
            c_v = k_BT_vals * dS_dkBT
        else:
            # Not enough points to compute derivative
            c_v = np.zeros_like(k_BT_vals)
            
        heat_cap_list[indices] = c_v

    df_reset["heat_capacity"] = heat_cap_list
    
    if is_multiindex:
        df_reset = df_reset.set_index(original_idx.names)
        
    df["heat_capacity"] = df_reset["heat_capacity"].values
    
    return df
