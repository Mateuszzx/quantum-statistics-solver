#include <cmath>
#include <stdexcept>
#include <algorithm>
#include "qss/math_utils.hpp"

namespace qss::math {


Mat make_empty_mat(std::size_t rows, std::size_t cols) {
  return Mat(rows, std::vector<double>(cols, 0.0));
}


// ---------- integrals ----------
double mean_abs_2d(const Mat& z, const Vec& kx, const Vec& ky) {
  const std::size_t ny = z.size();
  const std::size_t nx = (ny ? z[0].size() : 0);
  if (ny == 0 || nx == 0) return 0.0;

  // Mirror Python: use |z| and simple mean by number of samples
  long double sum = 0.0L;
  for (std::size_t j = 0; j < ny; ++j)
    for (std::size_t i = 0; i < nx; ++i)
      sum += std::fabsl(static_cast<long double>(z[j][i]));

  return static_cast<double>(sum / static_cast<long double>(kx.size() * ky.size()));
}

// ---------- bracketing: expansion ----------
bool find_bracket_expansion(const std::function<double(double)>& f,
                            double a0, double b0,
                            double expand_factor,
                            int max_expands,
                            double* a_out, double* b_out) {
  if (a0 >= b0) std::swap(a0, b0);
  double fa = f(a0);
  double fb = f(b0);
  if (std::isnan(fa) || std::isnan(fb)) return false;

  if (fa * fb < 0.0) {
    if (a_out) *a_out = a0;
    if (b_out) *b_out = b0;
    return true;
  }

  // Expand outward geometrically around the midpoint.
  const double mid = 0.5 * (a0 + b0);
  double radius = 0.5 * (b0 - a0);

  for (int k = 0; k < max_expands; ++k) {
    radius *= expand_factor;
    const double a = mid - radius;
    const double b = mid + radius;
    const double fa2 = f(a);
    const double fb2 = f(b);
    if (std::isnan(fa2) || std::isnan(fb2)) continue;
    if (fa2 * fb2 < 0.0) {
      if (a_out) *a_out = a;
      if (b_out) *b_out = b;
      return true;
    }
  }
  return false;
}



} // namespace qss::math
