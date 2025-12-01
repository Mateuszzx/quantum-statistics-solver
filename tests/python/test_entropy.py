
import numpy as np
import pandas as pd
import pytest
from unittest.mock import patch

from qss.entropy import entropy, entropy_single_channel, safe_log

def test_safe_log():
    arr = np.array([0.0, 1.0, np.e, -1.0])
    res = safe_log(arr)
    assert res[0] == 0.0
    assert res[1] == 0.0
    assert res[2] == pytest.approx(1.0)
    assert res[3] == 0.0  # Negative values handled as 0.0 by where condition

def test_entropy_single_channel_fd():
    n = np.array([0.5])
    # S = -0.5 ln 0.5 - 0.5 ln 0.5 = - ln 0.5 = ln 2
    s = entropy_single_channel(n, "FD")
    assert s[0] == pytest.approx(np.log(2))

def test_entropy_single_channel_flb():
    n = np.array([0.5])
    # S = S_FD + n ln 2 = ln 2 + 0.5 ln 2 = 1.5 ln 2
    s = entropy_single_channel(n, "FLB")
    assert s[0] == pytest.approx(1.5 * np.log(2))

def test_entropy_single_channel_nflb():
    # NFLB valid for n <= 0.5
    n = np.array([0.1])
    # S = -n ln n + n ln 2 + (1-n) ln(1-n) - (1-2n) ln(1-2n)
    s = entropy_single_channel(n, "NFLB")
    expected = -0.1 * np.log(0.1) + 0.1 * np.log(2) + 0.9 * np.log(0.9) - 0.8 * np.log(0.8)
    assert s[0] == pytest.approx(expected)

def test_entropy_single_channel_invalid():
    with pytest.raises(ValueError):
        entropy_single_channel(np.array([0.5]), "UNKNOWN")

@pytest.fixture
def sample_dataframe():
    index = pd.MultiIndex.from_tuples(
        [("FD", 0.5, 1.0, 0.0), ("FD", 0.5, 1.0, 1.0)],
        names=["distro", "n", "k_BT", "H_a"]
    )
    df = pd.DataFrame({"mu": [0.1, 0.2]}, index=index)
    return df

@pytest.fixture
def sample_energy_matrix():
    return np.array([[1.0, 2.0], [3.0, 4.0]])

def test_entropy_calculation_ha_zero(sample_dataframe, sample_energy_matrix):
    # H_a = 0
    # Mock distribution to return 1.0 (n_tot = 1.0 => n_single = 0.5)
    # S_single(0.5) = ln 2
    # Total S = 2 * ln 2
    with patch("qss.entropy.distribution", return_value=0.5 * np.ones_like(sample_energy_matrix)):
        df_result = entropy(sample_dataframe.iloc[[0]], sample_energy_matrix)
        assert "entropy" in df_result.columns
        assert df_result["entropy"].iloc[0] == pytest.approx(2 * np.log(2))

def test_entropy_calculation_ha_nonzero(sample_dataframe, sample_energy_matrix):
    # H_a != 0
    # Mock distribution to return 1.0 (n_tot = 1.0 => n_up = 0.5, n_down = 0.5)
    # S_up = ln 2, S_down = ln 2
    # Total S = 2 * ln 2
    with patch("qss.entropy.distribution", return_value=0.5 * np.ones_like(sample_energy_matrix)):
        df_result = entropy(sample_dataframe.iloc[[1]], sample_energy_matrix)
        assert "entropy" in df_result.columns
        assert df_result["entropy"].iloc[0] == pytest.approx(2 * np.log(2))

def test_entropy_flat_columns(sample_energy_matrix):
    df = pd.DataFrame({
        "distro": ["FD"],
        "n": [0.5],
        "k_BT": [1.0],
        "H_a": [0.0],
        "mu": [0.1]
    })
    with patch("qss.entropy.distribution", return_value=0.5 * np.ones_like(sample_energy_matrix)):
        df_result = entropy(df, sample_energy_matrix)
        assert "entropy" in df_result.columns
        assert df_result["entropy"].iloc[0] == pytest.approx(2 * np.log(2))

def test_entropy_malformed_dataframe(sample_energy_matrix):
    # DataFrame missing required index/columns
    df = pd.DataFrame({"wrong_col": [1, 2, 3]})

    with pytest.raises(ValueError, match="Input DataFrame must have"):
        entropy(df, sample_energy_matrix)


def test_entropy_ssl_calculation(sample_energy_matrix):
    # Test SSL entropy calculation
    # n_up = 0.25, n_down = 0.25 => n_hole = 0.5
    # S = -0.25 ln 0.25 - 0.25 ln 0.25 - 0.5 ln 0.5

    df = pd.DataFrame({
        "distro": ["SSL"],
        "n": [0.5],
        "k_BT": [1.0],
        "H_a": [0.0],
        "mu": [0.1]
    })

    n_val = 0.25
    expected_S = -2 * n_val * np.log(n_val) - (1 - 2*n_val) * np.log(1 - 2*n_val)

    with patch("qss.entropy.distribution", return_value=n_val * np.ones_like(sample_energy_matrix)):
        df_result = entropy(df, sample_energy_matrix)
        assert df_result["entropy"].iloc[0] == pytest.approx(expected_S)
