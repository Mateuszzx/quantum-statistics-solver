#include <stdexcept>
#include "qss/distributions.hpp"
#include "qss/math_utils.hpp"
#include "qss/solver.hpp"
#include "qss/types.hpp"
#include "chem_helpers.hpp"


namespace qss {


double chemical_potential(
    DistType dist,
    double beta,
    double n_assumed,
    const Vec& kx, const Vec& ky,
    const Mat& E,
    double delta,                     // = zeeman_coeff * B
    // bracketing/solver
    const IBracketer& bracketer, const BracketerOptions& bopt,
    const IRootSolver& solver,   const RootOptions& ropt,
    double a0 = -5000.0, double b0 = 1000.0)
{
  if (E.empty() || E[0].empty())
    throw std::invalid_argument("chemical_potential_2d_split: empty E");
  if (ky.size() != E.size() || kx.size() != E[0].size())
    throw std::invalid_argument("chemical_potential_2d_split: kx/ky vs E shape mismatch");

  auto objective = [&](double mu) -> double {
    Mat F = qss::chem::distribution_total(E, mu, beta, dist, delta);
    return qss::math::mean_abs_2d(F, kx, ky) - n_assumed;
  };

  Bracket br = bracketer.bracket(objective, a0, b0, bopt);
  if (!br.ok) throw std::runtime_error("failed to bracket chemical potential");
  return solver.solve(objective, br.a, br.b, ropt);
}

} // namespace qss


