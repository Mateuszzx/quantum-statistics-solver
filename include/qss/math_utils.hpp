#pragma once
#include <vector>
#include "types.hpp"

namespace qss::math {

// numerically-stable 1/(1+exp(x))
double logistic_neg(double x);

// map with saturation outside [down, up]
template <class F>
void saturating_map(const std::vector<double>& x,
                    double up_limit, double down_limit,
                    double value_when_big_pos,
                    double value_when_big_neg,
                    F&& f,
                    std::vector<double>& out);

// small helpers
Mat  make_empty_mat(std::size_t rows, std::size_t cols);

} // namespace qss::math

// ---- template impl (header-only) ----
namespace qss::math {
template <class F>
inline void saturating_map(const std::vector<double>& x,
                           double up_limit, double down_limit,
                           double value_when_big_pos,
                           double value_when_big_neg,
                           F&& f,
                           std::vector<double>& out)
{
  out.resize(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double xi = x[i];
    if (xi > up_limit)       out[i] = value_when_big_pos;
    else if (xi < down_limit)out[i] = value_when_big_neg;
    else                     out[i] = f(xi);
  }
}
} // namespace qss::math
