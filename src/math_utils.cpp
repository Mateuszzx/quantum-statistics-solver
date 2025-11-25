#include "qss/math_utils.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace qss::math
{


double mean_abs_2d(const Mat &z)
{
    const std::size_t ny = z.size();
    const std::size_t nx = (ny ? z[0].size() : 0);
    if (ny == 0 || nx == 0)
        return 0.0;

    long double sum = 0.0L;
    for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t i = 0; i < nx; ++i)
            sum += std::fabs(static_cast<long double>(z[j][i]));

    // Use actual sample count (ny*nx). The passed N may be lattice size but
    // relying on ny*nx avoids mismatch / scaling instability.
    return static_cast<double>(sum / static_cast<long double>(ny * nx));
}

} // namespace qss::math
