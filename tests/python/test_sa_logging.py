import os
import tempfile
import pytest
from qss.solver_options import sa_defaults
from qss.chemical_potential import chemical_potential
from qss.distributions import DistType
import numpy as np

def test_sa_logging():
    with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
        log_file = tmp.name

    try:
        # Setup a simple problem
        # We'll use chemical_potential with SA solver
        # It doesn't matter if it converges well, we just want to check logging

        # Mock data
        N = 10
        E = np.zeros((N, N)) # Flat band

        opts = sa_defaults(
            max_iters=10,
            log_file=log_file
        )

        # Run solver
        # This might fail to converge or whatever, but it should run SA
        try:
            chemical_potential(
                distros=[DistType.FD],
                n_vals=[0.5],
                kBT_vals=[1.0],
                H_a_vals=[0.0],
                N=N,
                E=E,
                solver_opts=opts.to_qss()
            )
        except Exception as e:
            print(f"Solver failed with: {e}")
            # Ignore solver errors, we check log file
            pass

        # Check if log file exists and has content
        assert os.path.exists(log_file)
        with open(log_file, 'r') as f:
            lines = f.readlines()

        # Header + at least 1 iteration
        assert len(lines) >= 2
        assert "iteration,temperature,energy" in lines[0]

        # Check first data line
        data = lines[1].split(',')
        assert len(data) >= 4 # iter, T, E, x0...
        assert int(data[0]) == 1 # iteration 1

    finally:
        if os.path.exists(log_file):
            os.remove(log_file)
