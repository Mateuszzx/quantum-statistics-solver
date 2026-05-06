"""Small Streamlit helpers (progress, matplotlib figures)."""

from __future__ import annotations

from typing import Any, Callable

import matplotlib.pyplot as plt


def st_progress_pair(st_mod: Any) -> tuple[Callable[[int, int], None], Any]:
    """Return ``(callback(done,total), placeholder)`` for experiment progress."""
    bar = st_mod.progress(0, text="Simulations: 0 / ?")

    def _cb(done: int, total: int) -> None:
        if total <= 0:
            return
        bar.progress(min(1.0, done / total), text=f"Simulations: {done} / {total}")

    return _cb, bar


def show_current_figure(st_mod: Any, *, clear_before: bool = True) -> None:
    if clear_before:
        fig = plt.gcf()
        if fig.get_axes():
            st_mod.pyplot(fig, clear_figure=True)
        return
    st_mod.pyplot(plt.gcf(), clear_figure=True)

