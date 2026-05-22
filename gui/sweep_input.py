"""Streamlit widgets for flexible parameter sweeps (scalar / list / grid)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import streamlit as st
from numpy.typing import NDArray

from gui.sweep_grids import values_from_input

SWEEP_MODES = ("Scalar", "List", "Linspace", "Logspace")


@dataclass(frozen=True)
class SweepParamSpec:
    """Default values for one swept quantity in the Panels view."""

    label: str
    scalar_default: float
    list_default: str
    grid_min: float
    grid_max: float
    grid_points: int
    mode_default: str = "Linspace"
    scalar_min: float | None = None
    scalar_max: float | None = None
    scalar_step: float | None = None


def parameter_input(
    *,
    param_key: str,
    spec: SweepParamSpec,
    container: st.delta_generator.DeltaGenerator,
) -> NDArray[np.float64]:
    """Render sweep controls for one parameter; return the value array."""
    container.markdown(spec.label)
    mode = container.radio(
        "Input mode",
        SWEEP_MODES,
        index=SWEEP_MODES.index(spec.mode_default)
        if spec.mode_default in SWEEP_MODES
        else 0,
        horizontal=True,
        key=f"{param_key}_mode",
    )

    scalar = spec.scalar_default
    list_text = spec.list_default
    vmin, vmax, count = spec.grid_min, spec.grid_max, spec.grid_points

    if mode == "Scalar":
        kwargs: dict[str, object] = {"value": spec.scalar_default, "key": f"{param_key}_scalar"}
        if spec.scalar_min is not None:
            kwargs["min_value"] = spec.scalar_min
        if spec.scalar_max is not None:
            kwargs["max_value"] = spec.scalar_max
        if spec.scalar_step is not None:
            kwargs["step"] = spec.scalar_step
        scalar = float(container.number_input("Value", **kwargs))
    elif mode == "List":
        list_text = container.text_input(
            "Values (comma-separated)",
            value=spec.list_default,
            key=f"{param_key}_list",
        )
    else:
        c1, c2 = container.columns(2)
        vmin = float(c1.number_input("Min", value=spec.grid_min, key=f"{param_key}_min"))
        vmax = float(c2.number_input("Max", value=spec.grid_max, key=f"{param_key}_max"))
        count = int(
            container.number_input(
                "Points",
                min_value=1,
                max_value=500,
                value=spec.grid_points,
                step=1,
                key=f"{param_key}_count",
            )
        )
        if mode == "Logspace" and count >= 2 and (vmin <= 0 or vmax <= 0):
            container.error("Log spacing needs strictly positive min and max.")

    try:
        vals = values_from_input(
            mode,
            scalar=scalar,
            list_text=list_text,
            vmin=vmin,
            vmax=vmax,
            count=count,
        )
    except ValueError as e:
        container.error(str(e))
        vals = np.array([float(spec.scalar_default)], dtype=np.float64)

    preview = ", ".join(f"{v:.6g}" for v in vals[:8])
    if len(vals) > 8:
        preview += f", … ({len(vals)} total)"
    container.caption(f"**{len(vals)}** value(s): {preview}")
    return vals


PANEL_SWEEP_SPECS: dict[str, SweepParamSpec] = {
    "n": SweepParamSpec(
        label=r"**Band filling** $n$",
        scalar_default=0.5,
        list_default="0.1, 0.3, 0.5, 0.7, 0.9",
        grid_min=0.1,
        grid_max=0.99,
        grid_points=20,
        mode_default="Linspace",
        scalar_min=0.001,
        scalar_max=1.0,
        scalar_step=0.01,
    ),
    "k_BT": SweepParamSpec(
        label=r"**Temperature** $k_B T / |t|$",
        scalar_default=0.1,
        list_default="0.01, 0.05, 0.1, 0.5, 1.0",
        grid_min=0.01,
        grid_max=1.0,
        grid_points=4,
        mode_default="Logspace",
    ),
    "H_a": SweepParamSpec(
        label=r"**Field** $H_a$ (T)",
        scalar_default=0.0,
        list_default="0.0, 1.0, 2.0",
        grid_min=0.0,
        grid_max=0.0,
        grid_points=1,
        mode_default="Scalar",
    ),
    "U": SweepParamSpec(
        label=r"**Hubbard** $U / |t|$",
        scalar_default=0.0,
        list_default="0.0, 2.0, 4.0, 8.0",
        grid_min=0.0,
        grid_max=0.0,
        grid_points=1,
        mode_default="Scalar",
    ),
}

def experiment_parameter_sweeps(
    container: st.delta_generator.DeltaGenerator,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Configure n, k_BT, H_a, and U sweeps for the Panels experiment."""
    sweeps: dict[str, NDArray[np.float64]] = {}
    for key, spec in PANEL_SWEEP_SPECS.items():
        title = key if key != "k_BT" else r"$k_B T/|t|$"
        with container.expander(title, expanded=False) as exp:
            sweeps[key] = parameter_input(
                param_key=f"sweep_{key}",
                spec=spec,
                container=exp,
            )
    return sweeps["n"], sweeps["k_BT"], sweeps["H_a"], sweeps["U"]
