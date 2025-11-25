#include "bindings.hpp"
#include "qss/solver.hpp"
#include <nanobind/nanobind.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

namespace nb = nanobind;

namespace qss_bindings
{

void bind_solver(nb::module_ &m)
{
    // Enums
    nb::enum_<qss::SolverType>(m, "SolverType")
        .value("Bisection", qss::SolverType::Bisection)
        .value("Secant", qss::SolverType::Secant)
        .value("Brent", qss::SolverType::Brent)
        .value("SimulatedAnnealing", qss::SolverType::SimulatedAnnealing);

    nb::enum_<qss::BracketerType>(m, "BracketerType")
        .value("Expansion", qss::BracketerType::Expansion)
        .value("Bisection", qss::BracketerType::Bisection);

    nb::enum_<qss::CoolingSchedule>(m, "CoolingSchedule")
        .value("Exponential", qss::CoolingSchedule::Exponential)
        .value("Linear", qss::CoolingSchedule::Linear)
        .value("Logarithmic", qss::CoolingSchedule::Logarithmic);

    // Options
    nb::class_<qss::SolverOptions>(m, "SolverOptions")
        .def(nb::init<>())
        .def_rw("abs_tol", &qss::SolverOptions::abs_tol)
        .def_rw("rel_tol", &qss::SolverOptions::rel_tol)
        .def_rw("max_iters", &qss::SolverOptions::max_iters);

    nb::class_<qss::BracketerOptions>(m, "BracketerOptions")
        .def(nb::init<>())
        .def_rw("expand_factor", &qss::BracketerOptions::expand_factor)
        .def_rw("max_expands", &qss::BracketerOptions::max_expands)
        .def_rw("initial_step", &qss::BracketerOptions::initial_step);

    nb::class_<qss::SimulatedAnnealingOptions>(m, "SimulatedAnnealingOptions")
        .def(nb::init<>())
        .def_rw("T0", &qss::SimulatedAnnealingOptions::T0)
        .def_rw("T_min", &qss::SimulatedAnnealingOptions::T_min)
        .def_rw("alpha", &qss::SimulatedAnnealingOptions::alpha)
        .def_rw("sched", &qss::SimulatedAnnealingOptions::sched)
        .def_rw("max_iters", &qss::SimulatedAnnealingOptions::max_iters)
        .def_rw("min_iters", &qss::SimulatedAnnealingOptions::min_iters)
        .def_rw("stall_iters", &qss::SimulatedAnnealingOptions::stall_iters)
        .def_rw("step_sigma", &qss::SimulatedAnnealingOptions::step_sigma)
        .def_rw("step_decay", &qss::SimulatedAnnealingOptions::step_decay)
        .def_rw("scale_sigma_by_T", &qss::SimulatedAnnealingOptions::scale_sigma_by_T)
        .def_rw("bounds", &qss::SimulatedAnnealingOptions::bounds)
        .def_rw("reflect_bounds", &qss::SimulatedAnnealingOptions::reflect_bounds)
        .def_rw("seed", &qss::SimulatedAnnealingOptions::seed);

    nb::class_<qss::CombinedSolverOptions>(m, "CombinedSolverOptions")
        .def(nb::init<>())
        .def_rw("solver", &qss::CombinedSolverOptions::solver)
        .def_rw("bracketer_type", &qss::CombinedSolverOptions::bracketer_type)
        .def_rw("bracketer_opts", &qss::CombinedSolverOptions::bracketer_opts)
        .def_rw("root", &qss::CombinedSolverOptions::root)
        .def_rw("sa", &qss::CombinedSolverOptions::sa);
}

} // namespace qss_bindings

