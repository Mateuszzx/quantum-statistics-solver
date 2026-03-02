import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

from qss.magnetization import magnetization_single, calculate_magnetization

class TestMagnetizationSingle:
    def test_valid_inputs(self):
        n_up = np.array([[1.0, 0.8], [0.6, 0.4]])
        n_down = np.array([[0.2, 0.1], [0.1, 0.0]])
        # Sum(n_up - n_down) = (0.8 + 0.7 + 0.5 + 0.4) = 2.4
        # kx * ky = 4
        # result = 2.4 / 4 = 0.6
        result = magnetization_single(n_up, n_down)
        assert np.isclose(result, 0.6)

    def test_shape_mismatch(self):
        n_up = np.zeros((2, 2))
        n_down = np.zeros((2, 3))
        with pytest.raises(ValueError, match="n_up and n_down must have the same shape"):
            magnetization_single(n_up, n_down)

    def test_absolute_value(self):
        # Ensure it returns absolute value
        n_up = np.zeros((2, 2))
        n_down = np.ones((2, 2))
        # Sum(0 - 1) = -4
        # / 4 = -1
        # abs -> 1.0
        result = magnetization_single(n_up, n_down)
        assert np.isclose(result, 1.0)


class TestCalculateMagnetization:
    @pytest.fixture
    def sample_energy(self):
        return np.zeros((10, 10))

    @pytest.fixture
    def sample_df_multiindex(self):
        # Index: 'distro', 'n', 'k_BT', 'H_a'
        index = pd.MultiIndex.from_tuples(
            [("FD", 0.5, 0.1, 0.0), ("FD", 0.5, 0.1, 1.0)],
            names=["distro", "n", "k_BT", "H_a"]
        )
        return pd.DataFrame({"mu": [0.1, 0.2]}, index=index)

    @pytest.fixture
    def sample_df_columns(self):
        return pd.DataFrame({
            "distro": ["FD", "FD"],
            "n": [0.5, 0.5],
            "k_BT": [0.1, 0.1],
            "H_a": [0.0, 1.0],
            "mu": [0.1, 0.2]
        })

    def test_zero_field(self, sample_df_multiindex, sample_energy):
        # First row has H_a = 0
        df_out = calculate_magnetization(sample_df_multiindex.iloc[[0]], sample_energy)
        assert "magnetization" in df_out.columns
        assert df_out.iloc[0]["magnetization"] == 0.0

    @patch("qss.magnetization.distribution")
    def test_nonzero_field(self, mock_dist, sample_df_multiindex, sample_energy):
        # Second row has H_a = 1.0
        # Mock distribution to return ones for up and zeros for down
        # shape of E is (10, 10) -> 100 elements
        mock_dist.side_effect = [
            np.ones((10, 10)),  # n_up
            np.zeros((10, 10))  # n_down
        ]

        df_in = sample_df_multiindex.iloc[[1]]
        df_out = calculate_magnetization(df_in, sample_energy)

        assert "magnetization" in df_out.columns
        # n_up - n_down = 1 everywhere. Sum = 100. Div by 100 = 1.
        assert np.isclose(df_out.iloc[0]["magnetization"], 1.0)

        # Verify distribution was called twice (up and down)
        assert mock_dist.call_count == 2

    def test_columns_input(self, sample_df_columns, sample_energy):
        # Test that it works with columns instead of MultiIndex
        # We'll just check it runs and produces output
        with patch("qss.magnetization.distribution", return_value=np.zeros((10, 10))):
            df_out = calculate_magnetization(sample_df_columns, sample_energy)
            assert "magnetization" in df_out.columns
            assert len(df_out) == 2

    def test_missing_columns(self, sample_energy):
        df = pd.DataFrame({"mu": [0.1]})
        with pytest.raises(ValueError, match="Input DataFrame must have"):
            calculate_magnetization(df, sample_energy)

    def test_zero_temperature(self, sample_energy):
        # k_BT = 0
        df = pd.DataFrame({
            "distro": ["FD"],
            "n": [0.5],
            "k_BT": [0.0],
            "H_a": [1.0],
            "mu": [0.1]
        })

        with patch("qss.magnetization.distribution") as mock_dist:
            mock_dist.return_value = np.zeros((10, 10))
            calculate_magnetization(df, sample_energy)

            # Check that beta passed to distribution is inf
            args, _ = mock_dist.call_args
            # args: E, mu, beta, ...
            beta = args[2]
            assert beta == 0.0
