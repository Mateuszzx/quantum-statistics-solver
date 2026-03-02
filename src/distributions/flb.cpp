/**
 * @file flb.cpp
 * @brief Implementation of the Fermi-Liquid-Like (Spin Liquid) distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

namespace
{
/**
 * @brief Core calculation for Fermi-Liquid-Like distribution.
 *
 * Computes \f$ f(x) = \frac{2}{e^x + 2} \f$ where \f$ x = \beta(E - \mu) \f$.
 * This distribution arises in certain spin liquid models.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 2].
 */
inline double flb_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 1.0;

    const double ex = 2.0 * std::exp(-x);
    return ex / (1 + ex);
}
} // namespace

double flb_distribution(double e_k, double mu, double kBT, double delta, int spin_shift) noexcept
{
    return flb_from_x((e_k + (spin_shift * delta) - mu) / kBT);
}


} // namespace qss::detail
