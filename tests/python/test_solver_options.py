import pytest
from unittest.mock import patch, MagicMock

try:
    from qss import _qss  # native module
except ImportError:
    pytest.skip("qss._qss native module not available", allow_module_level=True)

from qss.solver_options import (
    SolverType,
    BracketerType,
    CoolingSchedule,
    SolverOptions,
    BracketerOptions,
    SimulatedAnnealingOptions,
    CombinedSolverOptions,
    brent_defaults,
    sa_defaults,
)


def test_solver_options_to_from_qss():
    py = SolverOptions(abs_tol=1e-6, rel_tol=1e-4, max_iters=1234)
    nb = py.to_qss()
    assert isinstance(nb, _qss.SolverOptions)
    assert nb.abs_tol == pytest.approx(1e-6)
    assert nb.rel_tol == pytest.approx(1e-4)
    assert nb.max_iters == 1234
    back = SolverOptions.from_qss(nb)
    assert back == py


def test_bracketer_options_to_from_qss():
    py = BracketerOptions(expand_factor=3.0, max_expands=42, initial_step=0.05)
    nb = py.to_qss()
    assert nb.expand_factor == pytest.approx(3.0)
    assert nb.max_expands == 42
    assert nb.initial_step == pytest.approx(0.05)
    back = BracketerOptions.from_qss(nb)
    assert back == py


def test_sa_options_seed_none_roundtrip():
    py = SimulatedAnnealingOptions(
        T0=1.2,
        T_min=1e-7,
        alpha=0.9,
        sched=CoolingSchedule.Linear,
        max_iters=5555,
        min_iters=10,
        stall_iters=100,
        step_sigma=0.15,
        step_decay=0.98,
        scale_sigma_by_T=False,
        bounds=[(-5.0, 5.0)],
        reflect_bounds=False,
        seed=None,
    )
    nb = py.to_qss()
    assert isinstance(nb, _qss.SimulatedAnnealingOptions)
    assert nb.T0 == pytest.approx(1.2)
    assert nb.sched == _qss.CoolingSchedule.Linear
    back = SimulatedAnnealingOptions.from_qss(nb)
    # The backend assigns a random seed if None is provided, so we expect an integer back
    assert isinstance(back.seed, int)
    assert back.T0 == pytest.approx(py.T0)


def test_sa_options_seed_set_roundtrip():
    py = SimulatedAnnealingOptions(seed=12345)
    nb = py.to_qss()
    assert int(getattr(nb, "seed")) == 12345
    back = SimulatedAnnealingOptions.from_qss(nb)
    assert back.seed == 12345


def test_combined_brent_roundtrip():
    py = brent_defaults(
        abs_tol=1e-6,
        rel_tol=1e-4,
        max_iters=6000,
        bracketer=BracketerType.Expansion,
        expand_factor=3.0,
        max_expands=96,
        initial_step=0.05,
    )
    nb = py.to_qss()
    assert nb.solver == _qss.SolverType.Brent
    assert nb.bracketer_type == _qss.BracketerType.Expansion
    back = CombinedSolverOptions.from_qss(nb)
    assert back.solver == py.solver
    assert back.root.abs_tol == pytest.approx(1e-6)


def test_combined_sa_roundtrip():
    py = sa_defaults(T0=2.0, seed=999, bounds=[(-8.0, 8.0)], step_sigma=0.3)
    nb = py.to_qss()
    assert nb.solver == _qss.SolverType.SimulatedAnnealing
    assert nb.sa.T0 == pytest.approx(2.0)
    assert list(nb.sa.bounds) == [(-8.0, 8.0)]
    assert int(getattr(nb.sa, "seed")) == 999
    back = CombinedSolverOptions.from_qss(nb)
    assert back.sa.seed == 999


def test_solver_options_validation():
    """Test validation logic in SolverOptions."""
    with pytest.raises(ValueError, match="abs_tol must be positive"):
        SolverOptions(abs_tol=-1.0)
    with pytest.raises(ValueError, match="rel_tol must be positive"):
        SolverOptions(rel_tol=0.0)
    with pytest.raises(ValueError, match="max_iters must be positive"):
        SolverOptions(max_iters=0)


def test_bracketer_options_validation():
    """Test validation logic in BracketerOptions."""
    with pytest.raises(ValueError, match="expand_factor must be > 1.0"):
        BracketerOptions(expand_factor=1.0)
    with pytest.raises(ValueError, match="max_expands must be positive"):
        BracketerOptions(max_expands=0)
    with pytest.raises(ValueError, match="initial_step must be positive"):
        BracketerOptions(initial_step=-0.1)


def test_sa_options_validation():
    """Test validation logic in SimulatedAnnealingOptions."""
    with pytest.raises(ValueError, match="T0 must be positive"):
        SimulatedAnnealingOptions(T0=0.0)
    with pytest.raises(ValueError, match="T_min must be positive"):
        SimulatedAnnealingOptions(T_min=-1e-5)
    with pytest.raises(ValueError, match="T_min must be < T0"):
        SimulatedAnnealingOptions(T0=1.0, T_min=2.0)
    with pytest.raises(ValueError, match="alpha must be in"):
        SimulatedAnnealingOptions(alpha=1.1)
    with pytest.raises(ValueError, match="step_decay must be in"):
        SimulatedAnnealingOptions(step_decay=0.0)
    with pytest.raises(ValueError, match="bounds must be non-empty"):
        SimulatedAnnealingOptions(bounds=[])
    with pytest.raises(ValueError, match="Invalid bounds"):
        SimulatedAnnealingOptions(bounds=[(5.0, 2.0)])


def test_cooling_schedule_enums():
    """Test all CoolingSchedule enum values."""
    for sched in CoolingSchedule:
        nb = sched.to_qss()
        assert isinstance(nb, _qss.CoolingSchedule)
        back = CoolingSchedule.from_qss(nb)
        assert back == sched


def test_solver_type_enums():
    """Test all SolverType enum values."""
    for solver in SolverType:
        nb = solver.to_qss()
        assert isinstance(nb, _qss.SolverType)
        back = SolverType.from_qss(nb)
        assert back == solver


def test_sa_options_from_qss_zero_seed():
    """Test that seed=0 or -1 from C++ is converted to None."""
    # Mock a native object
    nb = MagicMock(spec=_qss.SimulatedAnnealingOptions)
    nb.T0 = 1.0
    nb.T_min = 1e-6
    nb.alpha = 0.95
    nb.sched = _qss.CoolingSchedule.Exponential
    nb.max_iters = 100
    nb.min_iters = 10
    nb.stall_iters = 20
    nb.step_sigma = 0.1
    nb.step_decay = 0.99
    nb.scale_sigma_by_T = True
    nb.bounds = [(-1.0, 1.0)]
    nb.reflect_bounds = True

    # Case 1: seed = 0
    nb.seed = 0
    py = SimulatedAnnealingOptions.from_qss(nb)
    assert py.seed is None

    # Case 2: seed = -1
    nb.seed = -1
    py = SimulatedAnnealingOptions.from_qss(nb)
    assert py.seed is None


def test_to_qss_runtime_error_when_module_missing():
    """Test that to_qss raises RuntimeError if _qss is None."""
    # We need to patch the module-level _qss variable in solver_options
    with patch("qss.solver_options._qss", None):
        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            SolverOptions().to_qss()

        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            BracketerOptions().to_qss()

        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            SimulatedAnnealingOptions().to_qss()

        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            CombinedSolverOptions().to_qss()

        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            SolverType.Brent.to_qss()

        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            BracketerType.Expansion.to_qss()

        with pytest.raises(RuntimeError, match="C\\+\\+ extension not available"):
            CoolingSchedule.Linear.to_qss()
