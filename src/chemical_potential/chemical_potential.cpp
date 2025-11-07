// src/chemical_potential.cpp
#include "qss/chemical_potential.hpp"
#include "chem_helpers.hpp" // distribution_total(E, mu, beta, dist, delta)
#include "qss/math_utils.hpp"
#include "qss/solver.hpp"
#include <stdexcept>

namespace qss
{

namespace
{
inline double zeeman_delta(double H_a_T, double g, double muB_eV_per_T)
{
    return 0.5 * g * muB_eV_per_T * H_a_T; // same units as E
}

static double solve_one_mu(DistType dist,
                           double beta,
                           double n_assumed,
                           const double N,
                           const Mat &E,
                           double delta,
                           const CombinedSolverOptions &opts,
                           double a0 = -5.0,
                           double b0 = 5.0,
                           double x0 = 0.0)
{
    if (E.empty() || E[0].empty())
        throw std::invalid_argument("chemical_potential_B: empty E");
    if (N != E.size() || N != E[0].size())
        throw std::invalid_argument(
            "chemical_potential_B: kx/ky vs E shape mismatch");

    // Objective with sign (for root solvers). |.| handled by solve_scalar when SA is used.
    auto f = [&](double mu) -> double {
        Mat F = qss::chem::distribution_total(E, mu, beta, dist, delta);
        return qss::math::mean_abs_2d(F, N) - n_assumed;
    };

    qss::OptimizationSolver solver(opts);
    return solver.solve_scalar(f, a0, b0, x0);
}
} // namespace

double chemical_potential_B(DistType dist,
                            double beta,
                            double n_assumed,
                            const double N,
                            const Mat &E,
                            double H_a_T,
                            double g,
                            double muB_eV_per_T,
                            const CombinedSolverOptions &solver_opts,
                            double a0,
                            double b0,
                            double x0)
{
    const double delta = zeeman_delta(H_a_T, g, muB_eV_per_T);
    return solve_one_mu(dist,
                        beta,
                        n_assumed,
                        N,
                        E,
                        delta,
                        solver_opts,
                        a0,
                        b0,
                        x0);
}

Vec chemical_potential_B(DistType dist,
                         const Vec &betas,
                         double n_assumed,
                         const double N,
                         const Mat &E,
                         double H_a_T,
                         double g,
                         double muB_eV_per_T,
                         const CombinedSolverOptions &solver_opts,
                         double a0,
                         double b0,
                         double x0)
{
    const double delta = zeeman_delta(H_a_T, g, muB_eV_per_T);
    Vec mus;
    mus.reserve(betas.size());
    for (double b : betas)
        mus.push_back(solve_one_mu(dist,
                                   b,
                                   n_assumed,
                                   N,
                                   E,
                                   delta,
                                   solver_opts,
                                   a0,
                                   b0,
                                   x0));
    return mus;
}

} // namespace qss
