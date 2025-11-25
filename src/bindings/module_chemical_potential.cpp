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
