/**
 * @file nfl.cpp
 * @brief Implementation of the Non-Fermi-Liquid distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

namespace
{
/**1
 * @brief Core calculation for Non-Fermi-Liquid distribution.
 *
 * Computes \f$ f(x) = \frac{1}{1 + 0.5e^x + \sqrt{1 + (e^x / 2 )^2}} \f$ where \f$ x = \beta(E - \mu) \f$.
 * This distribution arises in certain spin liquid models.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 1].
 */
inline double nfl_from_x(double x, double x_prim) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 0.5;

    const double ex = std::exp(-x);
    const double ex_prim = std::exp(-x_prim);

    return 0.5 + ex - 0.5 * std::sqrt(1.0 + 4.0 * ex_prim * ex_prim);
}
} // namespace

double nfl_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept
{
    return nfl_from_x(beta * (e_k + (spin_shift * delta) - mu), beta * (e_k -(spin_shift * delta) - mu));
}


} // namespace qss::detail
