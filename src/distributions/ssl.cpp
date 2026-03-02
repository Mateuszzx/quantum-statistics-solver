/**
 * @file ssl.cpp
 * @brief Implementation of ... distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{


namespace
{
/**
 * @brief Core calculation for .... distribution.
 *
 * Computes \f$ f(x) = \frac{e^x}{e^x + e^x'prim + 1} \f$ where \f$ x = \beta(E - \mu) \f$.
 * Handles numerical stability for large positive and negative x.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 0.5].
 */
inline double ssl_from_x(double x, double x_prim) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 1. / (1. + std::exp(x-x_prim));

    const double ex = std::exp(-x);
    const double ex_prim = std::exp(-x_prim);

    return (ex)/ (1.0 + ex + ex_prim);
}
} // namespace

double ssl_distribution(double e_k, double mu, double kBT, double delta, int spin_shift) noexcept
{
    return ssl_from_x((e_k + (spin_shift * delta) - mu) / kBT, (e_k - (spin_shift * delta) - mu) / kBT);
}

} // namespace qss::detail
