import numpy as np
import pytest

from qss.chemical_potential import chemical_potential
from qss.distributions import DistType, prepare_2D_lattice_data
from qss.solver_options import brent_defaults


def test_fd_mu_with_correct_t_prime_dispersion_reference_case():
    """Regression test for the diagonal t_prime tight-binding dispersion."""
    N = 100
    t = -1.0
    t_prim = 0.25
    kx = np.linspace(-np.pi, np.pi, N, endpoint=False)
    ky = np.linspace(-np.pi, np.pi, N, endpoint=False)

    e_k, _, _ = prepare_2D_lattice_data(kx, ky, t=t, t_prim=t_prim)

    mu_df = chemical_potential(
        distros=[DistType.FD],
        n_vals=[0.8],
        kBT_vals=[0.01],
        H_a_vals=[0.0],
        N=N,
        E=e_k,
        solver_opts=brent_defaults(),
        x0=-1.0,
        a0=-5.0,
        b0=5.0,
    )

    row = mu_df.reset_index().iloc[0]
    assert row["converged"]
    assert row["mu"] == pytest.approx(-0.9761154651641846, abs=1e-6)
