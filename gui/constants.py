"""Shared GUI constants (distributions, lattice defaults)."""

from __future__ import annotations

from qss.distributions import DistType

DIM = "2D"
DIM_2D = 2

DIST_MAP: dict[str, DistType] = {
    "FD": DistType.FD,
    "FLB": DistType.FLB,
    "NFL": DistType.NFL,
    "NFLB": DistType.NFLB,
    "SSL": DistType.SSL,
    "HUBBARD": DistType.HUBBARD,
}

SWEEP_VAR_LABELS: dict[str, str] = {
    "n": r"$n$",
    "mu": r"$\mu/|t|$",
    "k_BT": r"$k_B T/|t|$",
    "H_a": r"$H_a$",
    "U": r"$U/|t|$",
}
