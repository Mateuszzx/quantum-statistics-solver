"""Registry of experiment-derived quantities for the Panels GUI.

To add a quantity computed from ``mu_df`` and dispersion ``E``:

1. Implement the function in ``qss/`` (e.g. ``qss.my_quantity``).
2. Register it here with :func:`register_derived_quantity` (one call).
3. No changes needed in ``web_app.py``.

Multiple plot options can share one DataFrame (e.g. energy and free energy) by
reusing ``session_key`` and setting ``compute`` only on the first entry.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
import pandas as pd

from qss.energy import energy as compute_energy
from qss.entropy import entropy as compute_entropy
from qss.heat_capacity import heat_capacity as compute_heat_capacity
from qss.magnetization import calculate_magnetization
from qss.nkd import n_kd as compute_nkd

ComputeFn = Callable[[pd.DataFrame, np.ndarray], pd.DataFrame]


@dataclass(frozen=True)
class DerivedQuantitySpec:
    """One selectable quantity in sweep panels."""

    id: str
    display_name: str
    column: str
    ylabel: str
    session_key: str
    compute: ComputeFn | None = None


_REGISTRY: list[DerivedQuantitySpec] = []


def register_derived_quantity(spec: DerivedQuantitySpec) -> DerivedQuantitySpec:
    """Register a quantity; return *spec* for convenient chaining."""
    if any(q.id == spec.id for q in _REGISTRY):
        msg = f"Duplicate derived quantity id: {spec.id!r}"
        raise ValueError(msg)
    _REGISTRY.append(spec)
    return spec


def derived_quantities() -> tuple[DerivedQuantitySpec, ...]:
    """All registered quantities in registration order."""
    return tuple(_REGISTRY)


def derived_display_names() -> list[str]:
    return [q.display_name for q in _REGISTRY]


def derived_by_display_name(display_name: str) -> DerivedQuantitySpec:
    for q in _REGISTRY:
        if q.display_name == display_name:
            return q
    msg = f"Unknown quantity: {display_name!r}"
    raise KeyError(msg)


def compute_derived_dataframes(
    mu_df: pd.DataFrame,
    e_k: np.ndarray,
) -> dict[str, pd.DataFrame]:
    """Run each registered ``compute`` once per ``session_key``."""
    out: dict[str, pd.DataFrame] = {}
    for spec in _REGISTRY:
        if spec.compute is None or spec.session_key in out:
            continue
        out[spec.session_key] = spec.compute(mu_df, e_k)
    return out


def store_derived_in_session(
    session_state: object,
    *,
    mu_df: pd.DataFrame,
    e_k: np.ndarray,
) -> None:
    """Write ``mu_df`` and all post-processed DataFrames into Streamlit session."""
    session_state["mu_df"] = mu_df
    for key, df in compute_derived_dataframes(mu_df, e_k).items():
        session_state[key] = df


def derived_quantities_ready(session_state: object) -> bool:
    """True when every registered post-process quantity has been computed."""
    for spec in _REGISTRY:
        if spec.compute is not None and spec.session_key not in session_state:
            return False
    return True


def df_for_quantity(
    session_state: object,
    *,
    display_name: str,
) -> tuple[pd.DataFrame, str, str]:
    """Return *(dataframe, y_column, y_label)* for panel plotting."""
    spec = derived_by_display_name(display_name)
    df = session_state[spec.session_key]
    return df, spec.column, spec.ylabel


def derived_compute_summary() -> str:
    """Short label for progress spinners."""
    names = [q.display_name for q in _REGISTRY if q.compute is not None]
    return ", ".join(names)


# --- registrations (add new quantities below) ---

register_derived_quantity(
    DerivedQuantitySpec(
        id="mu",
        display_name="Chemical potential μ",
        column="mu",
        ylabel=r"chemical potential $\mu / |t|$",
        session_key="mu_df",
    )
)
register_derived_quantity(
    DerivedQuantitySpec(
        id="entropy",
        display_name="Entropy S",
        column="entropy",
        ylabel=r"entropy $S$",
        session_key="entropy_df",
        compute=compute_entropy,
    )
)
register_derived_quantity(
    DerivedQuantitySpec(
        id="heat_capacity",
        display_name="Heat capacity C",
        column="heat_capacity",
        ylabel=r"heat capacity $C$",
        session_key="heat_cap_df",
        compute=compute_heat_capacity,
    )
)
register_derived_quantity(
    DerivedQuantitySpec(
        id="energy",
        display_name="Energy E",
        column="energy",
        ylabel=r"energy $E$",
        session_key="energy_df",
        compute=compute_energy,
    )
)
register_derived_quantity(
    DerivedQuantitySpec(
        id="free_energy",
        display_name="Free energy F",
        column="free_energy",
        ylabel=r"free energy $F$",
        session_key="energy_df",
    )
)
register_derived_quantity(
    DerivedQuantitySpec(
        id="magnetization",
        display_name="Magnetization M",
        column="magnetization",
        ylabel=r"magnetization $M$",
        session_key="mag_df",
        compute=calculate_magnetization,
    )
)
register_derived_quantity(
    DerivedQuantitySpec(
        id="nkd",
        display_name="n_kd",
        column="nkd",
        ylabel=r"$n_{kd}$",
        session_key="nkd_df",
        compute=compute_nkd,
    )
)
