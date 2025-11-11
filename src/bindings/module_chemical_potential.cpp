#include <cstring>
#include <stdexcept>
#include <limits>
#include <vector>
#include <tuple>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/pair.h>

#include "bindings.hpp"
#include "qss/chemical_potential.hpp"
#include "qss/distributions.hpp"
#include "qss/solver_params.hpp"
#include "qss/types.hpp"

namespace nb = nanobind;
using nb::arg;

// Enforce C-contiguous arrays to ensure raw pointer access is safe
using Mat2 = nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>;

// Helper to safely copy numpy matrix to qss::Mat using raw pointers
static qss::Mat numpy_to_mat_safe(Mat2 E, size_t nx) {
    qss::Mat Em(nx, qss::Vec(nx));
    const double* ptr = E.data();
    for (size_t i = 0; i < nx; ++i) {
        for (size_t j = 0; j < nx; ++j) {
            Em[i][j] = ptr[i * nx + j];
        }
    }
    return Em;
}

namespace qss_bindings {

void bind_chemical_potential(nb::module_ &m) {
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

    nb::class_<qss::ChemicalPotentialResult>(m, "ChemicalPotentialResult")
        .def(nb::init<>())
        .def_rw("mu", &qss::ChemicalPotentialResult::mu)
        .def_rw("error", &qss::ChemicalPotentialResult::error)
        .def_rw("iterations", &qss::ChemicalPotentialResult::iterations)
        .def_rw("converged", &qss::ChemicalPotentialResult::converged);

    // Scalar binding ONLY - No GIL release, simple and safe
    m.def(
        "chemical_potential_B_scalar",
        [](qss::DistType dist, double n_assumed, double beta, double H_a_T, int N, Mat2 E,
           double g, double muB_eV_per_T, const qss::CombinedSolverOptions &opts,
           double a0, double b0, double x0) {

            size_t nx = static_cast<size_t>(N);
            if (E.shape(0) != nx || E.shape(1) != nx)
                throw std::invalid_argument("E shape mismatch");

            // Convert matrix (Safe copy)
            qss::Mat Em = numpy_to_mat_safe(E, nx);

            // Call C++ logic directly
            auto R = qss::chemical_potential_B(dist, beta, n_assumed, static_cast<double>(N),
                                               Em, H_a_T, g, muB_eV_per_T, opts, a0, b0, x0);

            return nb::make_tuple(R.mu, R.iterations, R.converged, R.error);
        },
        arg("distro"), arg("n_assumed"), arg("beta"), arg("H_a_T"), arg("N"), arg("E"),
        arg("g") = 2.0, arg("muB_eV_per_T") = qss::QSS_MU_B_EV_PER_T,
        arg("opts") = qss::CombinedSolverOptions(), arg("a0") = -10.0, arg("b0") = 10.0, arg("x0") = 0.0);

    m.def(
        "chemical_potential_B_multi_n",
        [](qss::DistType dist,
           double beta,
           const std::vector<double> &n_vals,
           double N,
           Mat2 E,
           double H_a,
           double g,
           double muB_eV_per_T,
           const qss::CombinedSolverOptions &opts,
           double a0,
           double b0,
           double x0) {
            auto Em = numpy_to_mat_safe(E, static_cast<size_t>(N));
            const auto results = qss::chemical_potential_B_multi_n(
                dist,
                beta,
                n_vals,
                N,
                Em,
                H_a,
                g,
                muB_eV_per_T,
                opts,
                a0,
                b0,
                x0);
            nb::list out;
            for (const auto &R : results) {
                out.append(nb::make_tuple(R.mu, R.iterations, R.converged, R.error));
            }
            return out;
        },
        arg("dist"),
        arg("beta"),
        arg("n_vals"),
        arg("N"),
        arg("E"),
        arg("H_a"),
        arg("g"),
        arg("muB_eV_per_T"),
        arg("solver_opts"),
        arg("a0") = -10.0,
        arg("b0") = 10.0,
        arg("x0") = 0.0);
}

} // namespace qss_bindings
