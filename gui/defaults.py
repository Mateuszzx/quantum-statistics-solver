"""Load GUI defaults from ``defaults.yaml`` (override path via ``QSS_GUI_DEFAULTS``)."""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

_BUILTIN = Path(__file__).resolve().parent / "defaults.yaml"


@lru_cache(maxsize=1)
def gui_defaults() -> dict[str, Any]:
    """Parsed defaults dict; reload process to pick up file edits."""
    path = Path(os.environ.get("QSS_GUI_DEFAULTS", _BUILTIN))
    with path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if not isinstance(data, dict):
        msg = f"GUI defaults must be a mapping, got {type(data).__name__}"
        raise TypeError(msg)
    return data


def cfg(*keys: str) -> Any:
    """Navigate nested keys, e.g. ``cfg("lattice", "t")``."""
    node: Any = gui_defaults()
    for key in keys:
        if not isinstance(node, dict) or key not in node:
            msg = f"Missing GUI default: {'.'.join(keys)} (failed at {key!r})"
            raise KeyError(msg)
        node = node[key]
    return node


def option_index(options: list[str], value: str, *, field: str) -> int:
    """Index of *value* in *options*; raises if missing."""
    try:
        return options.index(value)
    except ValueError as e:
        msg = f"Invalid {field} in defaults.yaml: {value!r} (choices: {options})"
        raise ValueError(msg) from e
