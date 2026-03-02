
import numpy as np
import pytest
from unittest.mock import patch, MagicMock

try:
    from qss import _qss
    HAS_QSS = True
except ImportError:
    HAS_QSS = False

from qss.distributions import distribution, DistType, _get_cpp_dist_type

@pytest.mark.skipif(not HAS_QSS, reason="_qss module not available")
class TestDistributionWrapper:
    def test_get_cpp_dist_type(self):
        assert _get_cpp_dist_type("FD") == _qss.DistType.FD
        assert _get_cpp_dist_type("HUBBARD") == _qss.DistType.HUBBARD
        assert _get_cpp_dist_type(DistType.FD) == _qss.DistType.FD
        assert _get_cpp_dist_type(_qss.DistType.FD) == _qss.DistType.FD

        with pytest.raises(ValueError):
            _get_cpp_dist_type("INVALID")

    def test_distribution_1d(self):
        e_k = np.linspace(-1, 1, 10)
        # We can't easily check the result values without knowing the C++ implementation details,
        # but we can check it runs and returns correct shape.
        res = distribution(e_k, 0.0, 1.0, 0.0, 0, dist_type="FD")
        assert res.shape == e_k.shape
        assert np.all(res >= 0)
        assert np.all(res <= 2) # FD max occupancy is 2 (spin degenerate)

    def test_distribution_2d(self):
        e_k = np.random.rand(5, 5)
        res = distribution(e_k, 0.0, 1.0, 0.0, 0, dist_type=DistType.FD)
        assert res.shape == e_k.shape

    def test_distribution_invalid_shape(self):
        e_k = np.zeros((2, 2, 2))
        with pytest.raises(ValueError, match="must be 1D or 2D"):
            distribution(e_k, 0.0, 1.0, 0.0, 0, dist_type="FD")

    def test_distribution_types(self):
        e_k = np.array([0.0])
        for dt in ["FD", "FLB", "SSL", "NFL", "NFLB", "HUBBARD"]:
            res = distribution(e_k, 0.0, 1.0, 0.0, 0, dist_type=dt)
            assert res.shape == (1,)
