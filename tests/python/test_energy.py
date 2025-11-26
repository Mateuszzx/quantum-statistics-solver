
import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch, MagicMock

from qss.energy import energy

@pytest.fixture
def sample_dataframe():
    # Create a sample DataFrame with MultiIndex
    index = pd.MultiIndex.from_tuples(
        [("FD", 0.5, 1.0, 0.0), ("FD", 0.5, 1.0, 1.0)],
        names=["distro", "n", "k_BT", "H_a"]
    )
    df = pd.DataFrame({"mu": [0.1, 0.2], "entropy": [0.5, 0.6]}, index=index)
    return df

@pytest.fixture
def sample_dataframe_no_entropy():
    # Create a sample DataFrame without entropy column
    index = pd.MultiIndex.from_tuples(
        [("FD", 0.5, 1.0, 0.0)],
        names=["distro", "n", "k_BT", "H_a"]
    )
    df = pd.DataFrame({"mu": [0.1]}, index=index)
    return df

@pytest.fixture
def sample_energy_matrix():
    return np.array([[1.0, 2.0], [3.0, 4.0]])

def test_energy_calculation_ha_zero(sample_dataframe, sample_energy_matrix):
    # Test case where H_a = 0
    # Mock distribution to return all ones
    with patch("qss.energy.distribution", return_value=np.ones_like(sample_energy_matrix)) as mock_dist:
        df_result = energy(sample_dataframe.iloc[[0]], sample_energy_matrix)
        
        assert "energy" in df_result.columns
        assert "free_energy" in df_result.columns
        
        # If n=1 everywhere, avg energy = mean(E * 1) = mean([1,2,3,4]) = 2.5
        expected_energy = 2.5
        assert df_result["energy"].iloc[0] == pytest.approx(expected_energy)
        
        # Free energy = E - T*S = 2.5 - 1.0 * 0.5 = 2.0
        assert df_result["free_energy"].iloc[0] == pytest.approx(2.0)

def test_energy_calculation_ha_nonzero(sample_dataframe, sample_energy_matrix):
    # Test case where H_a != 0
    # Mock distribution to return all ones (so n_up = 0.5, n_down = 0.5)
    with patch("qss.energy.distribution", return_value=np.ones_like(sample_energy_matrix)) as mock_dist:
        df_result = energy(sample_dataframe.iloc[[1]], sample_energy_matrix)
        
        # H_a = 1.0
        # delta = 0.5 * g * muB * H_a
        # E_up = E - delta
        # E_down = E + delta
        # n_up = 1/2, n_down = 1/2
        # avg = mean( (E-delta)*0.5 + (E+delta)*0.5 ) = mean( 0.5E - 0.5delta + 0.5E + 0.5delta ) = mean(E)
        
        expected_energy = 2.5
        assert df_result["energy"].iloc[0] == pytest.approx(expected_energy)

def test_energy_calculates_entropy_if_missing(sample_dataframe_no_entropy, sample_energy_matrix):
    with patch("qss.energy.distribution", return_value=np.ones_like(sample_energy_matrix)), \
         patch("qss.energy.calculate_entropy") as mock_entropy:
        
        # Mock entropy to return df with entropy column
        def side_effect(df, *args, **kwargs):
            df_out = df.copy()
            df_out["entropy"] = [0.7]
            return df_out
        
        mock_entropy.side_effect = side_effect
        
        df_result = energy(sample_dataframe_no_entropy, sample_energy_matrix)
        
        mock_entropy.assert_called_once()
        assert "entropy" in df_result.columns
        assert df_result["entropy"].iloc[0] == 0.7
        assert "energy" in df_result.columns

def test_energy_flat_columns(sample_energy_matrix):
    # Test with flat columns instead of MultiIndex
    df = pd.DataFrame({
        "distro": ["FD"],
        "n": [0.5],
        "k_BT": [1.0],
        "H_a": [0.0],
        "mu": [0.1],
        "entropy": [0.5]
    })
    
    with patch("qss.energy.distribution", return_value=np.ones_like(sample_energy_matrix)):
        df_result = energy(df, sample_energy_matrix)
        assert "energy" in df_result.columns
        assert df_result["energy"].iloc[0] == pytest.approx(2.5)

def test_energy_malformed_dataframe(sample_energy_matrix):
    # DataFrame missing required index/columns but HAS entropy (so we skip entropy calc)
    df = pd.DataFrame({
        "wrong_col": [1, 2, 3],
        "entropy": [0.1, 0.2, 0.3]
    })
    
    with pytest.raises(ValueError, match="Input DataFrame must have"):
        energy(df, sample_energy_matrix)

