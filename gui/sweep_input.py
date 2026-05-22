"""Streamlit widgets for flexible parameter sweeps (scalar / list / grid)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import streamlit as st
from numpy.typing import NDArray

from gui.defaults import cfg
from gui.sweep_grids import values_from_input

SWEEP_MODES = ("Scalar", "List", "Linspace", "Logspace")

_SWEEP_LABELS: dict[str, str] = {
    "n": r"**Band filling** $n$",
    "k_BT": r"**Temperature** $k_B T / |t|$",
    "H_a": r"**Field** $H_a$ (T)",
    "U": r"**Hubbard** $U / |t|$",
}


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


def _sweep_spec_from_yaml(key: str) -> SweepParamSpec:
    s = cfg("sweeps", key)
    scalar_min = float(s["scalar_min"]) if "scalar_min" in s else None
    scalar_max = float(s["scalar_max"]) if "scalar_max" in s else None
    scalar_step = float(s["scalar_step"]) if "scalar_step" in s else None
    return SweepParamSpec(
        label=_SWEEP_LABELS[key],
        scalar_default=float(s["scalar"]),
        list_default=str(s["list"]),
        grid_min=float(s["grid_min"]),
        grid_max=float(s["grid_max"]),
        grid_points=int(s["grid_points"]),
        mode_default=str(s["mode"]),
        scalar_min=scalar_min,
        scalar_max=scalar_max,
        scalar_step=scalar_step,
    )


def panel_sweep_specs() -> dict[str, SweepParamSpec]:
    """Sweep specs built from ``defaults.yaml``."""
    return {key: _sweep_spec_from_yaml(key) for key in ("n", "k_BT", "H_a", "U")}


def parameter_input(
    *,
    param_key: str,
    spec: SweepParamSpec,
    container: st.delta_generator.DeltaGenerator,
    grid_points_max: int,
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
                max_value=int(grid_points_max),
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


def experiment_parameter_sweeps(
    container: st.delta_generator.DeltaGenerator,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Configure n, k_BT, H_a, and U sweeps for the Panels experiment."""
    specs = panel_sweep_specs()
    grid_max = int(cfg("sweeps", "grid_points_max"))
    sweeps: dict[str, NDArray[np.float64]] = {}
    for key, spec in specs.items():
        title = key if key != "k_BT" else r"$k_B T/|t|$"
        with container.expander(title, expanded=False) as exp:
            sweeps[key] = parameter_input(
                param_key=f"sweep_{key}",
                spec=spec,
                container=exp,
                grid_points_max=grid_max,
            )
    return sweeps["n"], sweeps["k_BT"], sweeps["H_a"], sweeps["U"]
