#include <cmath>
#include "qss/math_utils.hpp"

namespace qss::math {

double logistic_neg(double x) {
  if (x >= 0.0) return 1.0 / (1.0 + std::exp(x));
  const double ex = std::exp(-x);
  return ex / (1.0 + ex);
}

Mat make_empty_mat(std::size_t rows, std::size_t cols) {
  return Mat(rows, std::vector<double>(cols, 0.0));
}

} // namespace qss::math
