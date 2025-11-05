import numpy as np
import pandas as pd


def numpy4d_to_dataframe(mu: np.ndarray,
                         distros,
                         n: np.ndarray,
                         beta: np.ndarray,
                         H_a: np.ndarray,
                         distros_labels: list[str] | None = None) -> pd.DataFrame:
    """
    Convert mu array with shape (ND, NN, NB, NH) -> DataFrame with
    index = MultiIndex.from_product([distros, beta, H_a])
    columns = n
    """
    if mu.ndim != 4:
        raise ValueError("mu must be a 4-D array with shape (ND, NN, NB, NH)")

    ND, NN, NB, NH = mu.shape
    if NN != len(n) or ND != len(distros) or NB != len(beta) or NH != len(H_a):
        raise ValueError("Lengths of distros, n, beta, H_a must match mu.shape")

    # Optional: map DistType enums to readable labels
    if distros_labels is None:
        distros_labels = [d.name if hasattr(d, "name") else str(d) for d in distros]

    idx = pd.MultiIndex.from_product([distros_labels, list(beta), list(H_a)],
                                     names=["distributions", "beta", "H_a"])

    # Reorder axes to (ND, NB, NH, NN) then flatten first three into rows
    arr = mu.transpose(0, 2, 3, 1).reshape(-1, NN)

    # Ensure column names are e.g. rounded floats or strings as you prefer
    col_names = list(n.astype(float))
    df = pd.DataFrame(arr, index=idx, columns=col_names)
    return df
