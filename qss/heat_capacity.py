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
    r"""
    Compute specific heat capacity for each state in the DataFrame.
    C = T * (\partial S / \partial T) = k_B T * (\partial S / \partial k_B T)
    """
    if "entropy" not in df.columns:
        df = calculate_entropy(df, E, g=g, muB_eV_per_T=muB_eV_per_T)

    # 1. Create a working copy and save the exact original row order
    df_work = df.copy()
    df_work["_orig_order"] = np.arange(len(df_work))

    # 2. Flatten the index
    df_work = df_work.reset_index()

    if "U" not in df_work.columns:
        df_work["U"] = U

    # 3. Sort values. np.gradient strictly requires ascending coordinates!
    group_cols = ["distro", "n", "H_a", "U"]
    df_work = df_work.sort_values([*group_cols, "k_BT"])

    # 4. Helper function to compute the derivative and return a mapped Series
    def compute_cv(group: pd.DataFrame) -> pd.Series:
        k = group["k_BT"].values.astype(float)
        S = group["entropy"].values.astype(float)

        if len(k) > 1:
            # SAFETY CHECK: Prevent division by zero if duplicate k_BT exist
            # by adding a microscopic offset to duplicates
            if np.any(np.diff(k) == 0):
                k = k + np.arange(len(k)) * 1e-12

            cv = k * np.gradient(S, k)
        else:
            cv = np.zeros_like(k)

        # Return a Series tied specifically to this group's index
        return pd.Series(cv, index=group.index)

    # 5. Apply the calculation. Because it returns a Series, it slots perfectly into the column.
    df_work["heat_capacity"] = df_work.groupby(group_cols, group_keys=False).apply(compute_cv)

    # 6. Restore the exact original row order using our tracker
    df_work = df_work.sort_values("_orig_order").drop(columns=["_orig_order"])

    # 7. Restore the original index cleanly
    if isinstance(df.index, pd.MultiIndex):
        df_work = df_work.set_index(df.index.names)
    else:
        df_work.index = df.index

    return df_work
