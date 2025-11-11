#pragma once
#include <functional>
#include <vector>

#include "types.hpp"

namespace qss::math
{

// small helpers
Mat make_empty_mat(std::size_t rows, std::size_t cols);


// -------- Integrals (left-edge) --------
double mean_abs_2d(const Mat &z, const double N);

// -------- Bracketing: expansion + (optional) contraction --------
/** Try to find [a,b] with f(a)*f(b) < 0 by expanding outwards.
 *  Returns true and sets (a_out,b_out) on success.
 */
bool find_bracket_expansion(const std::function<double(double)> &f,
                            double a0,
                            double b0,
                            double expand_factor = 1.5,
                            int max_expands = 200,
                            double *a_out = nullptr,
                            double *b_out = nullptr);


} // namespace qss::math
