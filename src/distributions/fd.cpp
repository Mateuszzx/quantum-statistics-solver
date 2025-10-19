#include <cmath>
#include "qss/distributions.hpp"
#include "qss/math_utils.hpp"

namespace qss {

double fd_distribution(double e_k, double mu, double beta) {
  const double x = beta * (e_k - mu);
  if (x > UP_EXP_LIMIT)   return 0.0;
  if (x < DOWN_EXP_LIMIT) return 2.0;
  return 2.0 * math::logistic_neg(x); // spin factor 2
}

Vec fd_distribution(const Vec& e_k, double mu, double beta) {
  std::vector<double> x; x.reserve(e_k.size());
  for (double ek : e_k) x.push_back(beta * (ek - mu));
  Vec out;
  math::saturating_map(
    x, UP_EXP_LIMIT, DOWN_EXP_LIMIT,
    /*big+*/ 0.0, /*big-*/ 2.0,
    [](double xi){ return 2.0 * math::logistic_neg(xi); },
    out
  );
  return out;
}

} // namespace qss
