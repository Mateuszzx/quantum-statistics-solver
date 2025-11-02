#include <cmath>
#include "impl.hpp"
#include "qss/types.hpp"
#include "../config.hpp"   // QSS_OMP_THRESHOLD

namespace qss::detail {

// ---- Implicit NFLB_i ----
// Solve on f ∈ [0, 1]:
//   p0 f^3 + p1 f^2 + p2 f + p3 = 0
// with p0 = e^x, p1 = 2 e^x + 8, p2 = -16, p3 = 8
// Return 2 * f_root (so result is in [0, 2]).
// Saturation: x >> 0 -> 0, x << 0 -> 2.

namespace {
inline bool find_root_on_unit_interval(double p0, double p1, double p2, double p3,
                                       double& root) noexcept {
  auto poly = [&](double f) noexcept {
    return ((p0 * f + p1) * f + p2) * f + p3;
  };

  double a = 0.0, b = 1.0;
  double fa = poly(a), fb = poly(b);
  if (std::isnan(fa) || std::isnan(fb)) return false;
  if (fa == 0.0) { root = a; return true; }
  if (fb == 0.0) { root = b; return true; }
  if (fa * fb > 0.0) return false; // no sign change

  // Robust bisection
  for (int it = 0; it < 80; ++it) {
    const double m = 0.5 * (a + b);
    const double fm = poly(m);
    if (std::abs(fm) < 1e-12) { root = m; return true; }
    if (fa * fm <= 0.0) { b = m; fb = fm; }
    else { a = m; fa = fm; }
  }
  root = 0.5 * (a + b);
  return true;
}

inline double nflb_i_from_x(double x) noexcept {
  if (x > UP_EXP_LIMIT)   return 0.0; // -> 0
  if (x < DOWN_EXP_LIMIT) return 2.0; // -> 2 (spin-like scaling)

  const double ex = std::exp(x);
  const double p0 = ex;
  const double p1 = 2.0 * ex + 8.0;
  const double p2 = -16.0;
  const double p3 = 8.0;

  double froot = 0.0;
  if (find_root_on_unit_interval(p0, p1, p2, p3, froot) && froot >= 0.0 && froot <= 1.0) {
    return 2.0 * froot; // scale to [0, 2]
  }
  return 0.0;
}
} // namespace

double nflb_i_distribution(double e_k, double mu, double beta) noexcept {
  return nflb_i_from_x(beta * (e_k - mu));
}

void nflb_i_distribution_into(const double* ek, std::size_t size,
                              double mu, double beta, double* out) {
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static) if(size >= QSS_OMP_THRESHOLD)
#endif
  for (std::size_t i = 0; i < size; ++i) {
    out[i] = nflb_i_from_x(beta * (ek[i] - mu));
  }
}

} // namespace qss::detail
