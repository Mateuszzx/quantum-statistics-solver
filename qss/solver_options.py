"""Utilities to configure and translate solver options between Python and C++.

This module defines Python enums and dataclasses that mirror the native `_qss`
solver configuration. Each dataclass provides `to_qss()`/`from_qss()` helpers to
convert between Python and native nanobind types seamlessly.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum

from . import _qss

# ----------------------------
# Public Python-side Enums
# ----------------------------


class SolverType(Enum):
    """Root-finding/minimization algorithms available in the backend."""

    Bisection = "Bisection"
    Secant = "Secant"
    Brent = "Brent"
    SimulatedAnnealing = "SimulatedAnnealing"

    def to_qss(self) -> _qss.SolverType:
        """Convert this enum to the native `_qss.SolverType`."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        return getattr(_qss.SolverType, self.name)

    @classmethod
    def from_qss(cls, val: _qss.SolverType) -> SolverType:
        """Create a Python enum from native `_qss.SolverType`."""
        return cls[val.name]


class BracketerType(Enum):
    """Strategies for finding an initial bracket for a root."""

    Expansion = "Expansion"
    Bisection = "Bisection"

    def to_qss(self) -> _qss.BracketerType:
        """Convert to native `_qss.BracketerType`."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        return getattr(_qss.BracketerType, self.name)

    @classmethod
    def from_qss(cls, val: _qss.BracketerType) -> BracketerType:
        """Create from native `_qss.BracketerType`."""
        return cls[val.name]


class CoolingSchedule(Enum):
    """Temperature schedules for Simulated Annealing."""

    Exponential = "Exponential"
    Linear = "Linear"
    Logarithmic = "Logarithmic"

    def to_qss(self) -> _qss.CoolingSchedule:
        """Convert to native `_qss.CoolingSchedule`."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        return getattr(_qss.CoolingSchedule, self.name)

    @classmethod
    def from_qss(cls, val: _qss.CoolingSchedule) -> CoolingSchedule:
        """Create from native `_qss.CoolingSchedule`."""
        return cls[val.name]


# ----------------------------
# Public Python-side Dataclasses
# ----------------------------


@dataclass
class SolverOptions:
    """Generic root-solver tolerances and iteration limits.

    Parameters
    ----------
    abs_tol : float
        Absolute tolerance for convergence.
    rel_tol : float
        Relative tolerance for convergence.
    max_iters : int
        Maximum number of iterations.
    """

    abs_tol: float = 1e-8
    rel_tol: float = 1e-6
    max_iters: int = 5000

    def __post_init__(self) -> None:
        """Validate parameters after initialization."""
        if self.abs_tol <= 0:
            raise ValueError("abs_tol must be positive")
        if self.rel_tol <= 0:
            raise ValueError("rel_tol must be positive")
        if self.max_iters <= 0:
            raise ValueError("max_iters must be positive")

    def to_qss(self) -> _qss.SolverOptions:
        """Create a native `_qss.SolverOptions` instance."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        o = _qss.SolverOptions()
        o.abs_tol = float(self.abs_tol)
        o.rel_tol = float(self.rel_tol)
        o.max_iters = int(self.max_iters)
        return o

    @classmethod
    def from_qss(cls, o: _qss.SolverOptions) -> SolverOptions:
        """Build from a native `_qss.SolverOptions`."""
        return cls(abs_tol=o.abs_tol, rel_tol=o.rel_tol, max_iters=o.max_iters)


@dataclass
class BracketerOptions:
    """Options for the bracketing phase prior to root solving.

    Parameters
    ----------
    expand_factor : float
        Geometric growth factor when expanding a bracket.
    max_expands : int
        Maximum number of expansion steps.
    initial_step : float
        Initial step size used when seeding a bracket around a guess.
    """

    expand_factor: float = 2.0
    max_expands: int = 64
    initial_step: float = 0.1

    def __post_init__(self) -> None:
        """Validate parameters after initialization."""
        if self.expand_factor <= 1.0:
            raise ValueError("expand_factor must be > 1.0")
        if self.max_expands <= 0:
            raise ValueError("max_expands must be positive")
        if self.initial_step <= 0:
            raise ValueError("initial_step must be positive")

    def to_qss(self) -> _qss.BracketerOptions:
        """Create a native `_qss.BracketerOptions`."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        o = _qss.BracketerOptions()
        o.expand_factor = float(self.expand_factor)
        o.max_expands = int(self.max_expands)
        o.initial_step = float(self.initial_step)
        return o

    @classmethod
    def from_qss(cls, o: _qss.BracketerOptions) -> BracketerOptions:
        """Build from a native `_qss.BracketerOptions`."""
        return cls(
            expand_factor=o.expand_factor,
            max_expands=o.max_expands,
            initial_step=o.initial_step,
        )


@dataclass
class SimulatedAnnealingOptions:
    """Options specific to the Simulated Annealing (SA) optimizer.

    Parameters
    ----------
    T0 : float
        Initial temperature.
    T_min : float
        Minimum temperature; annealing stops when below this (after min_iters).
    alpha : float
        Exponential decay factor used when `sched=Exponential`.
    sched : CoolingSchedule
        Temperature schedule kind (Exponential | Linear | Logarithmic).
    max_iters, min_iters, stall_iters : int
        Global iteration cap, minimum iterations before early stop checks,
        and allowed iterations without improvement.
    step_sigma : float
        Base proposal step standard deviation.
    step_decay : float
        Per-iteration multiplicative decay applied cumulatively to step size.
    scale_sigma_by_T : bool
        If True, scales step size by sqrt(T/T0).
    bounds : list[tuple[float, float]]
        Per-dimension (lo, hi) bounds; +/-inf is allowed.
    reflect_bounds : bool
        If True, reflect proposals at bounds (else clamp).
    seed : int | None
        RNG seed; None lets backend pick a default.
    """

    T0: float = 1.0
    T_min: float = 1e-6
    alpha: float = 0.95
    sched: CoolingSchedule = CoolingSchedule.Exponential

    max_iters: int = 5000
    min_iters: int = 50
    stall_iters: int = 200

    step_sigma: float = 0.2
    step_decay: float = 0.99
    scale_sigma_by_T: bool = True  # noqa: N815

    bounds: list[tuple[float, float]] = field(default_factory=lambda: [(-10.0, 10.0)])
    reflect_bounds: bool = True
    seed: int | None = None

    def __post_init__(self) -> None:
        """Validate parameters after initialization."""
        if self.T0 <= 0:
            raise ValueError("T0 must be positive")
        if self.T_min <= 0:
            raise ValueError("T_min must be positive")
        if self.T_min >= self.T0:
            raise ValueError("T_min must be < T0")
        if not 0 < self.alpha < 1:
            raise ValueError("alpha must be in (0, 1)")
        if self.max_iters <= 0:
            raise ValueError("max_iters must be positive")
        if self.min_iters < 0:
            raise ValueError("min_iters must be non-negative")
        if self.stall_iters <= 0:
            raise ValueError("stall_iters must be positive")
        if self.step_sigma <= 0:
            raise ValueError("step_sigma must be positive")
        if not 0 < self.step_decay <= 1:
            raise ValueError("step_decay must be in (0, 1]")
        if not self.bounds:
            raise ValueError("bounds must be non-empty")
        for lo, hi in self.bounds:
            if lo >= hi:
                raise ValueError(f"Invalid bounds: ({lo}, {hi}), lo must be < hi")

    def to_qss(self) -> _qss.SimulatedAnnealingOptions:
        """Create a native `_qss.SimulatedAnnealingOptions`."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        o = _qss.SimulatedAnnealingOptions()
        o.T0 = float(self.T0)
        o.T_min = float(self.T_min)
        o.alpha = float(self.alpha)
        o.sched = self.sched.to_qss()
        o.max_iters = int(self.max_iters)
        o.min_iters = int(self.min_iters)
        o.stall_iters = int(self.stall_iters)
        o.step_sigma = float(self.step_sigma)
        o.step_decay = float(self.step_decay)
        o.scale_sigma_by_T = bool(self.scale_sigma_by_T)
        o.bounds = [(float(a), float(b)) for (a, b) in self.bounds]
        o.reflect_bounds = bool(self.reflect_bounds)
        if self.seed is not None:
            o.seed = int(self.seed)
        return o

    @classmethod
    def from_qss(cls, o: _qss.SimulatedAnnealingOptions) -> SimulatedAnnealingOptions:
        """Build from a native `_qss.SimulatedAnnealingOptions`."""
        seed_opt: int | None = None
        raw_seed = getattr(o, "seed", None)
        if isinstance(raw_seed, int) and raw_seed not in (0, -1):
            seed_opt = raw_seed
        return cls(
            T0=o.T0,
            T_min=o.T_min,
            alpha=o.alpha,
            sched=CoolingSchedule.from_qss(o.sched),
            max_iters=o.max_iters,
            min_iters=o.min_iters,
            stall_iters=o.stall_iters,
            step_sigma=o.step_sigma,
            step_decay=o.step_decay,
            scale_sigma_by_T=o.scale_sigma_by_T,
            bounds=list(o.bounds),
            reflect_bounds=o.reflect_bounds,
            seed=seed_opt,
        )


@dataclass
class CombinedSolverOptions:
    """Top-level configuration for the solver pipeline.

    Combines the bracketing method, root solver, and optional SA settings.
    """

    solver: SolverType = SolverType.Brent
    bracketer_type: BracketerType = BracketerType.Expansion

    bracketer_opts: BracketerOptions = field(default_factory=BracketerOptions)
    root: SolverOptions = field(default_factory=SolverOptions)
    sa: SimulatedAnnealingOptions = field(default_factory=SimulatedAnnealingOptions)

    def to_qss(self) -> _qss.CombinedSolverOptions:
        """Create a native `_qss.CombinedSolverOptions`."""
        if _qss is None:
            raise RuntimeError("C++ extension not available")
        o = _qss.CombinedSolverOptions()
        o.solver = self.solver.to_qss()
        o.bracketer_type = self.bracketer_type.to_qss()
        o.bracketer_opts = self.bracketer_opts.to_qss()
        o.root = self.root.to_qss()
        o.sa = self.sa.to_qss()
        return o

    @classmethod
    def from_qss(cls, o: _qss.CombinedSolverOptions) -> CombinedSolverOptions:
        """Build from a native `_qss.CombinedSolverOptions`."""
        return cls(
            solver=SolverType.from_qss(o.solver),
            bracketer_type=BracketerType.from_qss(o.bracketer_type),
            bracketer_opts=BracketerOptions.from_qss(o.bracketer_opts),
            root=SolverOptions.from_qss(o.root),
            sa=SimulatedAnnealingOptions.from_qss(o.sa),
        )


# ----------------------------
# Convenience builders
# ----------------------------


def brent_defaults(  # noqa: PLR0913
    abs_tol: float = 1e-8,
    rel_tol: float = 1e-6,
    max_iters: int = 5000,
    bracketer: BracketerType = BracketerType.Expansion,
    expand_factor: float = 2.0,
    max_expands: int = 64,
    initial_step: float = 0.1,
) -> CombinedSolverOptions:
    """Ready-to-use Brent configuration with sensible defaults.

    Parameters
    ----------
    abs_tol, rel_tol, max_iters : float, float, int
        Tolerances and iteration cap for Brent solver.
    bracketer : BracketerType
        Bracketing strategy (Expansion or Bisection).
    expand_factor, max_expands, initial_step : float, int, float
        Parameters controlling the bracketer behavior.

    Returns
    -------
    CombinedSolverOptions
        Fully populated options object for use with the backend.
    """
    return CombinedSolverOptions(
        solver=SolverType.Brent,
        bracketer_type=bracketer,
        bracketer_opts=BracketerOptions(
            expand_factor=expand_factor,
            max_expands=max_expands,
            initial_step=initial_step,
        ),
        sa=SimulatedAnnealingOptions(),
        root=SolverOptions(abs_tol=abs_tol, rel_tol=rel_tol, max_iters=max_iters),
    )


def sa_defaults(  # noqa: PLR0913
    T0: float = 1.0,
    T_min: float = 1e-6,
    alpha: float = 0.95,
    sched: CoolingSchedule = CoolingSchedule.Exponential,
    max_iters: int = 5000,
    min_iters: int = 50,
    stall_iters: int = 200,
    step_sigma: float = 0.2,
    step_decay: float = 0.99,
    scale_sigma_by_T: bool = True,
    bounds: list[tuple[float, float]] | None = None,
    reflect_bounds: bool = True,
    seed: int | None = None,
) -> CombinedSolverOptions:
    """Ready-to-use Simulated Annealing configuration with sensible defaults.

    Parameters
    ----------
    T0, T_min, alpha : float
        Temperature schedule parameters.
    sched : CoolingSchedule
        Kind of temperature schedule.
    max_iters, min_iters, stall_iters : int
        Global iteration cap, minimum iterations before early stop checks,
        and allowed iterations without improvement.
    step_sigma, step_decay : float
        Base proposal step size and per-iteration multiplicative decay.
    scale_sigma_by_T : bool
        If True, scales step size by sqrt(T/T0).
    bounds : list[(float, float)] | None
        Optional per-dimension bounds; defaults to [(-10, 10)].
    reflect_bounds : bool
        Reflect at bounds if True; clamp if False.
    seed : int | None
        RNG seed; None lets backend pick a default.

    Returns
    -------
    CombinedSolverOptions
        Fully populated options object configured for SA.
    """
    return CombinedSolverOptions(
        solver=SolverType.SimulatedAnnealing,
        sa=SimulatedAnnealingOptions(
            T0=T0,
            T_min=T_min,
            alpha=alpha,
            sched=sched,
            max_iters=max_iters,
            min_iters=min_iters,
            stall_iters=stall_iters,
            step_sigma=step_sigma,
            step_decay=step_decay,
            scale_sigma_by_T=scale_sigma_by_T,
            bounds=[(-10.0, 10.0)] if bounds is None else bounds,
            reflect_bounds=reflect_bounds,
            seed=seed,
        ),
    )


__all__ = [
    "BracketerOptions",
    "BracketerType",
    "CombinedSolverOptions",
    "CoolingSchedule",
    "SimulatedAnnealingOptions",
    "SolverOptions",
    "SolverType",
    "brent_defaults",
    "sa_defaults",
]
