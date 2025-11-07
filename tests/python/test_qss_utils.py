"""
Unit tests for qss Python utilities (data_format.py, prep.py)
"""

import numpy as np
import pytest

from qss.data_format import numpy4d_to_dataframe
from qss.prep import prepare_2D_lattice_data


class TestDataFormat:
    """Tests for qss.data_format module"""

    def test_numpy4d_to_dataframe_basic(self):
        """Test basic 4D array to DataFrame conversion"""
        # Create test data: (ND=2, NN=3, NB=2, NH=2)
        mu = np.random.rand(2, 3, 2, 2)
        distros = ["FermiDirac", "BoseEinstein"]
        n = np.array([0.5, 1.0, 1.5])
        beta = np.array([0.1, 0.2])
        H_a = np.array([0.0, 1.0])

        df = numpy4d_to_dataframe(mu, distros, n, beta, H_a)

        # Check shape
        assert df.shape == (8, 3)  # (2*2*2 = 8 rows, 3 columns)

        # Check index structure
        assert df.index.names == ["distributions", "beta", "H_a"]
        assert len(df.index.levels[0]) == 2  # 2 distributions

        # Check columns
        assert len(df.columns) == 3
        assert df.columns.tolist() == [0.5, 1.0, 1.5]

    def test_numpy4d_to_dataframe_with_labels(self):
        """Test conversion with custom distribution labels"""
        mu = np.ones((2, 2, 2, 2))

        class MockDistro:
            def __init__(self, name):
                self.name = name

        distros = [MockDistro("FD"), MockDistro("BE")]
        n = np.array([1.0, 2.0])
        beta = np.array([0.1, 0.2])
        H_a = np.array([0.0, 1.0])

        df = numpy4d_to_dataframe(mu, distros, n, beta, H_a)

        # Check that distribution names are extracted
        assert "FD" in df.index.get_level_values(0)
        assert "BE" in df.index.get_level_values(0)

    def test_numpy4d_to_dataframe_with_custom_labels(self):
        """Test conversion with explicit distribution labels"""
        mu = np.ones((2, 2, 2, 2))
        distros = [0, 1]  # Enum-like integers
        n = np.array([1.0, 2.0])
        beta = np.array([0.1, 0.2])
        H_a = np.array([0.0, 1.0])
        distros_labels = ["Fermi-Dirac", "Bose-Einstein"]

        df = numpy4d_to_dataframe(mu, distros, n, beta, H_a, distros_labels)

        # Check custom labels are used
        assert "Fermi-Dirac" in df.index.get_level_values(0)
        assert "Bose-Einstein" in df.index.get_level_values(0)

    def test_numpy4d_to_dataframe_wrong_dimensions(self):
        """Test that wrong dimensionality raises ValueError"""
        mu_3d = np.random.rand(2, 3, 2)  # Only 3D
        distros = ["FD", "BE"]
        n = np.array([0.5, 1.0, 1.5])
        beta = np.array([0.1, 0.2])
        H_a = np.array([0.0, 1.0])

        with pytest.raises(ValueError, match="mu must be a 4-D array"):
            numpy4d_to_dataframe(mu_3d, distros, n, beta, H_a)

    def test_numpy4d_to_dataframe_shape_mismatch(self):
        """Test that shape mismatch raises ValueError"""
        mu = np.random.rand(2, 3, 2, 2)
        distros = ["FD", "BE"]
        n = np.array([0.5, 1.0])  # Wrong length (should be 3)
        beta = np.array([0.1, 0.2])
        H_a = np.array([0.0, 1.0])

        with pytest.raises(ValueError, match="Lengths of distros, n, beta, H_a must match"):
            numpy4d_to_dataframe(mu, distros, n, beta, H_a)

    def test_numpy4d_to_dataframe_values(self):
        """Test that values are correctly transposed and reshaped"""
        # Create a simple test case with known values
        mu = np.arange(24).reshape(2, 3, 2, 2)  # (ND=2, NN=3, NB=2, NH=2)
        distros = ["D0", "D1"]
        n = np.array([1.0, 2.0, 3.0])
        beta = np.array([0.1, 0.2])
        H_a = np.array([0.0, 1.0])

        df = numpy4d_to_dataframe(mu, distros, n, beta, H_a)

        # Check first row corresponds to mu[0, :, 0, 0]
        first_row = df.iloc[0].values
        expected_first_row = mu[0, :, 0, 0]
        np.testing.assert_array_equal(first_row, expected_first_row)


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

    def test_lattice_to_dataframe_workflow(self):
        """Test complete workflow from lattice preparation to DataFrame"""
        # Prepare lattice data
        kx = np.linspace(-np.pi, np.pi, 3)
        ky = np.linspace(-np.pi, np.pi, 3)
        e_k, KX, KY = prepare_2D_lattice_data(kx, ky)

        # Simulate some mu data with same k-space shape
        # Shape: (ND=2 distros, NN=3 filling, NB=2 beta, NH=3*3 k-points)
        # For this test, flatten k-space as if it were H_a dimension
        mu_data = np.random.rand(2, 3, 2, 9)

        distros = ["FD", "BE"]
        n = np.array([0.5, 1.0, 1.5])
        beta = np.array([0.1, 0.2])
        H_a = e_k.flatten()  # Use flattened energies as "field" values

        df = numpy4d_to_dataframe(mu_data, distros, n, beta, H_a)

        # Verify integration
        assert df.shape[0] == 2 * 2 * 9  # ND * NB * NH
        assert df.shape[1] == 3  # NN
        assert df.index.names == ["distributions", "beta", "H_a"]


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


@pytest.mark.parametrize(
    "shape",
    [
        (2, 3, 2, 2),
        (1, 5, 3, 4),
        (3, 2, 4, 5),
    ],
)
def test_dataframe_with_different_shapes(shape):
    """Test DataFrame conversion with various array shapes"""
    ND, NN, NB, NH = shape
    mu = np.random.rand(*shape)

    distros = [f"D{i}" for i in range(ND)]
    n = np.linspace(0.1, 1.0, NN)
    beta = np.linspace(0.1, 1.0, NB)
    H_a = np.linspace(0.0, 10.0, NH)

    df = numpy4d_to_dataframe(mu, distros, n, beta, H_a)

    assert df.shape == (ND * NB * NH, NN)
    assert len(df.index.levels[0]) == ND
    assert len(df.columns) == NN
