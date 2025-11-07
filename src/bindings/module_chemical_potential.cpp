#include <cstring>
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/pair.h> // add: for vector<pair<double,double>>
#include <nanobind/stl/vector.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "bindings.hpp"
#include "qss/chemical_potential.hpp"
#include "qss/distributions.hpp"
#include "qss/solver.hpp"
#include "qss/types.hpp"

namespace nb = nanobind;
using nb::arg;

using Vec1 = nb::ndarray<nb::numpy, const double, nb::ndim<1>, nb::c_contig>;
using Mat2 = nb::ndarray<nb::numpy, const double, nb::ndim<2>, nb::c_contig>;
using Out4D = nb::ndarray<nb::numpy, double, nb::ndim<4>>;

static qss::SolverParams make_params(double tol,
                                     int max_iters,
                                     double expand_factor,
                                     int max_expands,
                                     double a0,
                                     double b0)
{
    qss::SolverParams p{};
    p.tol = tol;
    p.max_iters = max_iters;
    p.expand_factor = expand_factor;
    p.max_expands = max_expands;
    p.a0 = a0;
    p.b0 = b0;
    return p;
}

static Out4D make_out_4d(size_t d0, size_t d1, size_t d2, size_t d3)
{
    const size_t N = d0 * d1 * d2 * d3;
    double *data = new double[N];
    nb::capsule owner(data, [](void *p) noexcept { delete[] (double *)p; });
    return Out4D(data,
                 {d0, d1, d2, d3},
                 owner,
                 {(int64_t)(d1 * d2 * d3), (int64_t)(d2 * d3), (int64_t)d3, 1});
}

namespace qss_bindings
{

void bind_chemical_potential(nb::module_ &m)
{
    // Enums (unchanged)
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

    // Options: switch def_readwrite -> def_rw
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
        .def_rw("scale_sigma_by_T",
                &qss::SimulatedAnnealingOptions::scale_sigma_by_T)
        .def_rw(
            "bounds",
            &qss::SimulatedAnnealingOptions::bounds) // List[Tuple[float,float]]
        .def_rw("reflect_bounds",
                &qss::SimulatedAnnealingOptions::reflect_bounds)
        .def_rw("seed", &qss::SimulatedAnnealingOptions::seed);

    nb::class_<qss::CombinedSolverOptions>(m, "CombinedSolverOptions")
        .def(nb::init<>())
        .def_rw("solver", &qss::CombinedSolverOptions::solver)
        .def_rw("bracketer_type", &qss::CombinedSolverOptions::bracketer_type)
        .def_rw("bracketer_opts", &qss::CombinedSolverOptions::bracketer_opts)
        .def_rw("root", &qss::CombinedSolverOptions::root)
        .def_rw("sa", &qss::CombinedSolverOptions::sa);

    // ---------- grid API using CombinedSolverOptions ----------
    auto grid_lambda = [](std::vector<qss::DistType> distros,
                          Vec1 n_assumed,
                          Vec1 beta,
                          Vec1 H_a_T,
                          const int N,
                          Mat2 E,
                          double g,
                          double muB_eV_per_T,
                          const qss::CombinedSolverOptions &opts,
                          double a0,
                          double b0,
                          double x0,
                          int threads) {
        const size_t nx = (size_t)N;
        if ((size_t)E.shape(0) != nx || (size_t)E.shape(1) != nx)
            throw std::runtime_error("E must have shape (len(ky), len(kx))");

        const size_t ND = distros.size();
        const size_t NN = (size_t)n_assumed.shape(0);
        const size_t NB = (size_t)beta.shape(0);
        const size_t NH = (size_t)H_a_T.shape(0);
        if (!ND || !NN || !NB || !NH)
            throw std::runtime_error(
                "distros, n_assumed, beta, H_a_T must be non-empty");

        // Convert grid arrays once
        qss::Mat Em(nx, qss::Vec(nx));
        for (size_t j = 0; j < nx; ++j)
        {
            const double *row = E.data() + j * nx;
            std::memcpy(Em[j].data(), row, nx * sizeof(double));
        }

        auto out = make_out_4d(ND, NN, NB, NH);
        double *__restrict P = out.data();

        const size_t NBxNH = NB * NH;

        nb::gil_scoped_release release;
#ifdef _OPENMP
        if (threads > 0)
            omp_set_num_threads(threads);
#pragma omp parallel for collapse(3) schedule(static)
#endif
        for (ptrdiff_t id = 0; id < (ptrdiff_t)ND; ++id)
            for (ptrdiff_t inx = 0; inx < (ptrdiff_t)NN; ++inx)
                for (ptrdiff_t ib = 0; ib < (ptrdiff_t)NB; ++ib)
                {
                    const double beta_v = beta.data()[ib];
                    const double n_v = n_assumed.data()[inx];
                    const size_t base = (size_t)id * (NN * NBxNH) +
                                        (size_t)inx * NBxNH + (size_t)ib * NH;

                    for (ptrdiff_t ih = 0; ih < (ptrdiff_t)NH; ++ih)
                    {
                        const double mu =
                            qss::chemical_potential_B(distros[(size_t)id],
                                                      beta_v,
                                                      n_v,
                                                      nx,
                                                      Em,
                                                      H_a_T.data()[ih],
                                                      g,
                                                      muB_eV_per_T,
                                                      opts,
                                                      a0,
                                                      b0,
                                                      x0);
                        P[base + (size_t)ih] = mu;
                    }
                }
        return out;
    };

    // Keep API name users expect
    m.def("chemical_potential_B",
          grid_lambda,
          arg("distros"),
          arg("n_assumed"),
          arg("beta"),
          arg("H_a_T"),
          arg("N"),
          arg("E"),
          arg("g") = 2.0,
          arg("muB_eV_per_T") = qss::QSS_MU_B_EV_PER_T,
          arg("opts") = qss::CombinedSolverOptions(),
          arg("a0") = -10.0,
          arg("b0") = 10.0,
          arg("x0") = 0.0,
          arg("threads") = 0);

    // Backward-compatible alias (same impl)
    m.def("chemical_potential_B_grid",
          grid_lambda,
          arg("distros"),
          arg("n_assumed"),
          arg("beta"),
          arg("H_a_T"),
          arg("N"),
          arg("E"),
          arg("g") = 2.0,
          arg("muB_eV_per_T") = qss::QSS_MU_B_EV_PER_T,
          arg("opts") = qss::CombinedSolverOptions(),
          arg("a0") = -10.0,
          arg("b0") = 10.0,
          arg("x0") = 0.0,
          arg("threads") = 0);
}

} // namespace qss_bindings
