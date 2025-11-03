// src/chemical_potential.cpp
#include <stdexcept>
#include "qss/chemical_potential.hpp"
#include "qss/solver.hpp"
#include "qss/math_utils.hpp"
#include "chem_helpers.hpp"   // distribution_total(E, mu, beta, dist, delta)

namespace qss {

namespace {
inline double zeeman_delta(double H_a_T, double g, double muB_eV_per_T) {
  return 0.5 * g * muB_eV_per_T * H_a_T; // same units as E
}

static double solve_one_mu(
    DistType dist, double beta, double n_assumed,
    const Vec& kx, const Vec& ky, const Mat& E,
    double delta, const SolverParams& p)
{
  if (E.empty() || E[0].empty())
    throw std::invalid_argument("chemical_potential_B: empty E");
  if (ky.size() != E.size() || kx.size() != E[0].size())
    throw std::invalid_argument("chemical_potential_B: kx/ky vs E shape mismatch");

  auto br = make_bracketer(p.bracketer);
  auto sl = make_solver(p.solver);
  BracketerOptions bopt{p.expand_factor, p.max_expands};
  RootOptions      ropt{p.tol, p.max_iters};

  auto objective = [&](double mu) -> double {
    Mat F = qss::chem::distribution_total(E, mu, beta, dist, delta);
    return qss::math::mean_abs_2d(F, kx, ky) - n_assumed;
  };

  Bracket B = br->bracket(objective, p.a0, p.b0, bopt);
  if (!B.ok) throw std::runtime_error("chemical_potential_B: failed to bracket root");
  return sl->solve(objective, B.a, B.b, ropt);
}
} // namespace

double chemical_potential_B(
    DistType dist, double beta, double n_assumed,
    const Vec& kx, const Vec& ky, const Mat& E,
    double H_a_T, double g, double muB_eV_per_T,
    const SolverParams& params)
{
  const double delta = zeeman_delta(H_a_T, g, muB_eV_per_T);
  return solve_one_mu(dist, beta, n_assumed, kx, ky, E, delta, params);
}

Vec chemical_potential_B(
    DistType dist, const Vec& betas, double n_assumed,
    const Vec& kx, const Vec& ky, const Mat& E,
    double H_a_T, double g, double muB_eV_per_T,
    const SolverParams& params)
{
  const double delta = zeeman_delta(H_a_T, g, muB_eV_per_T);
  Vec mus; mus.reserve(betas.size());
  for (double b : betas)
    mus.push_back(solve_one_mu(dist, b, n_assumed, kx, ky, E, delta, params));
  return mus;
}

} // namespace qss
