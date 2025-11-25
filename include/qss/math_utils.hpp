#pragma once
#include <vector>

#include "types.hpp"

namespace qss::math
{

/**
 * @brief Calculation of 2D integral from left edges
 * 
 * @param z - matrix 
 * @return double 
 */
double mean_abs_2d(const Mat &z);

} // namespace qss::math
