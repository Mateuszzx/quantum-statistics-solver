# import pytest

# try:
#     from qss import _qss  # native module
# except Exception:
#     pytest.skip("qss._qss native module not available", allow_module_level=True)

# from qss.solver_options import (
#     SolverType,
#     BracketerType,
#     CoolingSchedule,
#     SolverOptions,
#     BracketerOptions,
#     SimulatedAnnealingOptions,
#     CombinedSolverOptions,
#     brent_defaults,
#     sa_defaults,
# )


# def test_solver_options_to_from_qss():
#     py = SolverOptions(abs_tol=1e-6, rel_tol=1e-4, max_iters=1234)
#     nb = py.to_qss()
#     assert isinstance(nb, _qss.SolverOptions)
#     assert nb.abs_tol == pytest.approx(1e-6)
#     assert nb.rel_tol == pytest.approx(1e-4)
#     assert nb.max_iters == 1234
#     back = SolverOptions.from_qss(nb)
#     assert back == py


# def test_bracketer_options_to_from_qss():
#     py = BracketerOptions(expand_factor=3.0, max_expands=42, initial_step=0.05)
#     nb = py.to_qss()
#     assert nb.expand_factor == pytest.approx(3.0)
#     assert nb.max_expands == 42
#     assert nb.initial_step == pytest.approx(0.05)
#     back = BracketerOptions.from_qss(nb)
#     assert back == py


# def test_sa_options_seed_none_roundtrip():
#     py = SimulatedAnnealingOptions(
#         T0=1.2,
#         T_min=1e-7,
#         alpha=0.9,
#         sched=CoolingSchedule.Linear,
#         max_iters=5555,
#         min_iters=10,
#         stall_iters=100,
#         step_sigma=0.15,
#         step_decay=0.98,
#         scale_sigma_by_T=False,
#         bounds=[(-5.0, 5.0)],
#         reflect_bounds=False,
#         seed=None,
#     )
#     nb = py.to_qss()
#     assert isinstance(nb, _qss.SimulatedAnnealingOptions)
#     assert nb.T0 == pytest.approx(1.2)
#     assert nb.sched == _qss.CoolingSchedule.Linear
#     back = SimulatedAnnealingOptions.from_qss(nb)
#     assert back.seed is None
#     assert back.T0 == pytest.approx(py.T0)


# def test_sa_options_seed_set_roundtrip():
#     py = SimulatedAnnealingOptions(seed=12345)
#     nb = py.to_qss()
#     assert int(getattr(nb, "seed")) == 12345
#     back = SimulatedAnnealingOptions.from_qss(nb)
#     assert back.seed == 12345


# def test_combined_brent_roundtrip():
#     py = brent_defaults(
#         abs_tol=1e-6,
#         rel_tol=1e-4,
#         max_iters=6000,
#         bracketer=BracketerType.Expansion,
#         expand_factor=3.0,
#         max_expands=96,
#         initial_step=0.05,
#     )
#     nb = py.to_qss()
#     assert nb.solver == _qss.SolverType.Brent
#     assert nb.bracketer_type == _qss.BracketerType.Expansion
#     back = CombinedSolverOptions.from_qss(nb)
#     assert back.solver == py.solver
#     assert back.root.abs_tol == pytest.approx(1e-6)


# def test_combined_sa_roundtrip():
#     py = sa_defaults(T0=2.0, seed=999, bounds=[(-8.0, 8.0)], step_sigma=0.3)
#     nb = py.to_qss()
#     assert nb.solver == _qss.SolverType.SimulatedAnnealing
#     assert nb.sa.T0 == pytest.approx(2.0)
#     assert list(nb.sa.bounds) == [(-8.0, 8.0)]
#     assert int(getattr(nb.sa, "seed")) == 999
#     back = CombinedSolverOptions.from_qss(nb)
#     assert back.sa.seed == 999
