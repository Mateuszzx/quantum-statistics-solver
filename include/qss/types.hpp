#pragma once
#include <vector>

namespace qss
{

using Vec = std::vector<double>;
using Mat = std::vector<std::vector<double>>;

inline constexpr double UP_EXP_LIMIT = 20.0;    // exp(+20) ~ 4.85e8
inline constexpr double DOWN_EXP_LIMIT = -16.0; // exp(-16) ~ 1.1e-7

} // namespace qss
