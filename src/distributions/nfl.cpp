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
// TODO: implement actual NFL distribution
/**
 * @brief Core calculation for Non-Fermi-Liquid distribution.
 *
 * Computes \f$ f(x) = \frac{2}{e^x + 2} \f$ where \f$ x = \beta(E - \mu) \f$.
 * This distribution arises in certain spin liquid models.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 1].
 */
inline double nfl_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 1.0;

    const double ex = 2.0 * std::exp(-x);
    return ex / (1 + ex);
}
} // namespace

double nfl_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept
{
    return nfl_from_x(beta * (e_k + (spin_shift * delta) - mu));
}


} // namespace qss::detail
