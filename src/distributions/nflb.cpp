#include <cmath>
#include "impl.hpp"
#include "qss/types.hpp"
#include "../config.hpp"   // QSS_OMP_THRESHOLD

namespace qss::detail {

// ---- Non-Fermi-Liquid (NFLB) ----
// For x = beta*(e_k - mu):
//   f(x) = (e^x - sqrt(e^{2x} + 8 e^x) + 8) / (e^x + 8)
// Saturation: x >> 0  -> 0,  x << 0 -> 1

namespace {
inline double nflb_from_x(double x) noexcept {
  if (x > UP_EXP_LIMIT)   return 0.0;
  if (x < DOWN_EXP_LIMIT) return 1.0;

  const double ex = std::exp(x);
  // Numerically OK up to x≈20 (your clamp); avoids overflow.
  const double num = ex - std::sqrt(ex * ex + 8.0 * ex) + 8.0;
  const double den = ex + 8.0;
  return num / den; // in [0, 1]
}
} // namespace

double nflb_distribution(double e_k, double mu, double beta) noexcept {
  return nflb_from_x(beta * (e_k - mu));
}

void nflb_distribution(const double* ek, std::size_t size,
                            double mu, double beta, double* out) {
#if defined(_OPENMP)
  #pragma omp parallel for schedule(static) if(size >= QSS_OMP_THRESHOLD)
#endif
  for (std::size_t i = 0; i < size; ++i) {
    out[i] = nflb_from_x(beta * (ek[i] - mu));
  }
}

} // namespace qss::detail
