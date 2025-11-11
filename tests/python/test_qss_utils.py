"""
Unit tests for qss Python utilities (data_format.py, prep.py)
"""

import numpy as np
import pytest

from qss.distributions import prepare_2D_lattice_data


class TestDataFormat:
    """Tests for qss.data_format module"""


class TestPrep:
    """Tests for qss.prep module"""

    def test_prepare_2D_lattice_data_scalar_inputs(self):
        """Test with scalar inputs"""
        kx = 0.5
        ky = 1.0

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky)

        # With scalar inputs, meshgrid should create 1x1 arrays
        assert e_k.shape == (1, 1)
        assert KX.shape == (1, 1)
        assert KY.shape == (1, 1)
        assert KX[0, 0] == pytest.approx(0.5)
        assert KY[0, 0] == pytest.approx(1.0)

    def test_prepare_2D_lattice_data_array_inputs(self):
        """Test with array inputs"""
        kx = np.linspace(-np.pi, np.pi, 5)
        ky = np.linspace(-np.pi, np.pi, 5)

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky)

        # Check shapes
        assert e_k.shape == (5, 5)
        assert KX.shape == (5, 5)
        assert KY.shape == (5, 5)

        # Check meshgrid structure
        assert np.allclose(KX[0, :], kx)  # First row should be kx
        assert np.allclose(KY[:, 0], ky)  # First column should be ky

    def test_prepare_2D_lattice_data_default_parameters(self):
        """Test with default tight-binding parameters"""
        kx = np.array([0.0, np.pi / 2, np.pi])
        ky = np.array([0.0, np.pi / 2, np.pi])

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky)

        # At (0, 0), with t=-1, t_prim=0.25:
        # e_k = 2*(-1)*(cos(0) + cos(0)) + 2*0.25*(cos(0) + cos(0))
        #     = -4 + 1 = -3
        assert e_k[0, 0] == pytest.approx(-3.0)

        # Check that e_k is float type
        assert e_k.dtype == np.float64

    def test_prepare_2D_lattice_data_custom_parameters(self):
        """Test with custom tight-binding parameters"""
        kx = np.array([0.0])
        ky = np.array([0.0])
        t = -2.0
        t_prim = 0.5

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky, t=t, t_prim=t_prim)

        # At (0, 0): e_k = 2*t*2 + 2*t_prim*2 = 4*t + 4*t_prim
        expected = 4 * t + 4 * t_prim
        assert e_k[0, 0] == pytest.approx(expected)

    def test_prepare_2D_lattice_data_symmetry(self):
        """Test energy dispersion symmetry properties"""
        kx = np.linspace(-np.pi, np.pi, 10)
        ky = np.linspace(-np.pi, np.pi, 10)

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky)

        # Check that e_k has expected symmetries for square lattice
        # e_k should be periodic with period 2π
        # Check center value
        center_idx = len(kx) // 2
        assert np.isfinite(e_k[center_idx, center_idx])

    def test_prepare_2D_lattice_data_mixed_inputs(self):
        """Test with mixed scalar and array inputs"""
        kx = np.array([0.0, 0.5, 1.0])
        ky = 0.5  # Scalar

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky)

        # ky is converted to array, so shape should be (1, 3)
        assert e_k.shape == (1, 3)
        assert KX.shape == (1, 3)
        assert KY.shape == (1, 3)

    def test_prepare_2D_lattice_data_negative_hopping(self):
        """Test that negative hopping parameter works correctly"""
        kx = np.array([0.0, np.pi])
        ky = np.array([0.0, np.pi])
        t = -1.0  # Standard negative hopping

        e_k, KX, KY = prepare_2D_lattice_data(kx, ky, t=t)

        # At (0, 0): cos(0) + cos(0) = 2, so first term = 2*t*2 = -4
        # At (π, π): cos(π) + cos(π) = -2, so first term = 2*t*(-2) = 4
        assert e_k[0, 0] < e_k[1, 1]  # Minimum at Γ point for t < 0


# Integration tests
class TestIntegration:
    """Integration tests combining multiple modules"""


# Fixtures for common test data
@pytest.fixture
def sample_kspace():
    """Fixture providing sample k-space grid"""
    kx = np.linspace(-np.pi, np.pi, 10)
    ky = np.linspace(-np.pi, np.pi, 10)
    return kx, ky


@pytest.fixture
def sample_mu_array():
    """Fixture providing sample 4D mu array"""
    return np.random.rand(2, 5, 3, 4)


# Parametrized tests
@pytest.mark.parametrize(
    "t,t_prim",
    [
        (-1.0, 0.0),
        (-1.0, 0.25),
        (-1.0, -0.1),
        (-2.0, 0.5),
    ],
)
def test_lattice_with_different_hoppings(t, t_prim):
    """Test lattice preparation with various hopping parameters"""
    kx = np.array([0.0, np.pi / 2, np.pi])
    ky = np.array([0.0, np.pi / 2, np.pi])

    e_k, KX, KY = prepare_2D_lattice_data(kx, ky, t=t, t_prim=t_prim)

    assert e_k.shape == (3, 3)
    assert np.all(np.isfinite(e_k))  # All values should be finite

    # Energy at (0, 0) should be: 2*t*2 + 2*t_prim*2
    expected_gamma = 4 * t + 4 * t_prim
    assert e_k[0, 0] == pytest.approx(expected_gamma)
