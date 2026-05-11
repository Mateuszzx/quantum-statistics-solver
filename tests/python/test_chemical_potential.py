
import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

try:
    from qss import _qss
    HAS_QSS = True
except ImportError:
    HAS_QSS = False

if HAS_QSS:
    from qss.chemical_potential import chemical_potential, _coerce_native_distro, _dist_label

@pytest.mark.skipif(not HAS_QSS, reason="_qss module not available")
class TestChemicalPotential:
    def test_coerce_native_distro(self):
        assert _coerce_native_distro("FD") == _qss.DistType.FD
        with pytest.raises(TypeError):
            _coerce_native_distro(123)

    def test_dist_label(self):
        assert _dist_label("FD") == "FD"
        assert _dist_label(_qss.DistType.FD) == "FD"

    @patch("qss.chemical_potential._qss.chemical_potential_B_multi_n")
    def test_chemical_potential_dataframe_structure(self, mock_calc):
        # Setup mock return
        # Returns list of tuples: (mu, iters, conv, err)
        # We pass 2 n_vals, so it returns 2 results per call
        mock_calc.return_value = [
            (0.1, 10, True, 1e-6),
            (0.2, 12, True, 1e-6)
        ]

        E = np.zeros((2, 2))
        n_vals = [0.4, 0.6]
        kBT_vals = [1.0]
        H_a_vals = [0.0]
        distros = ["FD"]

        df = chemical_potential(
            distros=distros,
            n_vals=n_vals,
            kBT_vals=kBT_vals,
            H_a_vals=H_a_vals,
            E=E
        )

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert list(df.columns) == ["mu", "iterations", "converged", "residual"]
        assert df.index.names == ["distro", "n", "k_BT", "H_a"]

        # Check index values
        assert df.index[0] == ("FD", 0.4, 1.0, 0.0)
        assert df.index[1] == ("FD", 0.6, 1.0, 0.0)

    def test_chemical_potential_invalid_E(self):
        E = np.zeros((2, 3)) # Not square
        with pytest.raises(ValueError, match="square 2D array"):
            chemical_potential("FD", [0.5], [1.0], [0.0], E)

    @patch("qss.chemical_potential._qss.chemical_potential_B_multi_n")
    def test_chemical_potential_multiple_params(self, mock_calc):
        # 2 distros * 2 betas * 1 H_a * 1 n = 4 calls to backend?
        # No, backend is vectorized over n.
        # Loops: distros -> betas -> H_a.
        # Inside: call backend (handles all n).

        mock_calc.return_value = [(0.0, 0, True, 0.0)] # 1 n_val

        E = np.zeros((2, 2))
        n_vals = [0.5]
        kBT_vals = [1.0, 2.0]
        H_a_vals = [0.0]
        distros = ["FD", "FLB"]

        df = chemical_potential(
            distros=distros,
            n_vals=n_vals,
            kBT_vals=kBT_vals,
            H_a_vals=H_a_vals,
            E=E
        )

        # 2 distros * 2 betas * 1 H_a * 1 n = 4 rows
        assert len(df) == 4
        assert mock_calc.call_count == 4

    @patch("qss.chemical_potential._qss.chemical_potential_B_multi_n")
    def test_chemical_potential_with_solver_options(self, mock_calc):
        mock_calc.return_value = [(0.0, 0, True, 0.0)]
        E = np.zeros((2, 2))

        # Create custom solver options
        from qss.solver_options import SolverOptions, CombinedSolverOptions
        opts = CombinedSolverOptions(root=SolverOptions(abs_tol=1e-9))

        chemical_potential(
            distros="FD",
            n_vals=[0.5],
            kBT_vals=[1.0],
            H_a_vals=[0.0],
            E=E,
            solver_opts=opts
        )

        # Verify that the options were passed to the backend
        # The 9th argument to chemical_potential_B_multi_n is solver_opts
        args, _ = mock_calc.call_args
        passed_opts = args[9]
        assert isinstance(passed_opts, _qss.CombinedSolverOptions)
        assert passed_opts.root.abs_tol == pytest.approx(1e-9)

    @patch("qss.chemical_potential._qss.chemical_potential_B_multi_n")
    def test_chemical_potential_single_distro_string(self, mock_calc):
        mock_calc.return_value = [(0.0, 0, True, 0.0)]
        E = np.zeros((2, 2))

        df = chemical_potential(
            distros="FD", # Single string instead of list
            n_vals=[0.5],
            kBT_vals=[1.0],
            H_a_vals=[0.0],
            E=E
        )

        assert len(df) == 1
        assert df.index[0][0] == "FD"

    @patch("qss.chemical_potential._qss.chemical_potential_B_multi_n")
    def test_chemical_potential_custom_N(self, mock_calc):
        mock_calc.return_value = [(0.0, 0, True, 0.0)]
        N = 100
        E = np.zeros((N, N))

        chemical_potential(
            distros="FD",
            n_vals=[0.5],
            kBT_vals=[1.0],
            H_a_vals=[0.0],
            E=E,
            N=N # Custom system size
        )

        # Verify N was passed correctly (4th argument)
        args, _ = mock_calc.call_args
        assert args[3] == N

    def test_chemical_potential_N_must_match_E_shape(self):
        E = np.zeros((2, 2))
        with pytest.raises(ValueError, match="does not match N=100"):
            chemical_potential(
                distros="FD",
                n_vals=[0.5],
                kBT_vals=[1.0],
                H_a_vals=[0.0],
                E=E,
                N=100,
            )

    def test_chemical_potential_invalid_distro(self):
        E = np.zeros((2, 2))
        with pytest.raises(TypeError):
            chemical_potential(
                distros=[123], # Invalid type
                n_vals=[0.5],
                kBT_vals=[1.0],
                H_a_vals=[0.0],
                E=E
            )

    def test_coerce_native_distro_variants(self):
        # Test string
        assert _coerce_native_distro("FD") == _qss.DistType.FD
        # Test PyDistType
        from qss.distributions import DistType as PyDistType
        assert _coerce_native_distro(PyDistType.FD) == _qss.DistType.FD
        # Test native
        assert _coerce_native_distro(_qss.DistType.FD) == _qss.DistType.FD

    @patch("qss.chemical_potential._qss.chemical_potential_B_multi_n")
    def test_chemical_potential_native_solver_opts(self, mock_calc):
        mock_calc.return_value = [(0.0, 0, True, 0.0)]
        E = np.zeros((2, 2))

        # Create native solver options directly
        native_opts = _qss.CombinedSolverOptions()
        native_opts.root.abs_tol = 1e-7

        chemical_potential(
            distros="FD",
            n_vals=[0.5],
            kBT_vals=[1.0],
            H_a_vals=[0.0],
            E=E,
            solver_opts=native_opts
        )

        args, _ = mock_calc.call_args
        passed_opts = args[9]
        assert passed_opts is native_opts
        assert passed_opts.root.abs_tol == pytest.approx(1e-7)
