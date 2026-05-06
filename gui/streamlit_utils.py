"""Small Streamlit helpers (progress, matplotlib figures)."""

from __future__ import annotations

from typing import Callable, Protocol

import matplotlib.pyplot as plt


class StreamlitLike(Protocol):
    def progress(self, value: float, *, text: str | None = None) -> object: ...

    def pyplot(self, fig: object, *, clear_figure: bool = True) -> None: ...


def st_progress_pair(st_mod: StreamlitLike) -> tuple[Callable[[int, int], None], object]:
    """Return ``(callback(done,total), placeholder)`` for experiment progress."""
    bar = st_mod.progress(0, text="Simulations: 0 / ?")

    def _cb(done: int, total: int) -> None:
        if total <= 0:
            return
        # The Streamlit progress handle supports calling `.progress(...)` again.
        bar.progress(min(1.0, done / total), text=f"Simulations: {done} / {total}")

    return _cb, bar


def show_current_figure(st_mod: StreamlitLike, *, clear_before: bool = True) -> None:
    if clear_before:
        fig = plt.gcf()
        if fig.get_axes():
            st_mod.pyplot(fig, clear_figure=True)
        return
    st_mod.pyplot(plt.gcf(), clear_figure=True)
