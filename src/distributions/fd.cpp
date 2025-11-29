/**
 * @file fd.cpp
 * @brief Implementation of the Fermi-Dirac distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

namespace
{
/**
 * @brief Core calculation for Fermi-Dirac distribution.
 *
 * Computes \f$ f(x) = \frac{2}{e^x + 1} \f$ where \f$ x = \beta(E - \mu) \f$.
 * Handles numerical stability for large positive and negative x.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 1].
 */
inline double fd_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 1.0;
    if (x >= 0.0)
        return 1.0 / (1.0 + std::exp(x));
    const double ex = std::exp(-x);
    return ex / (1.0 + ex);
}
} // namespace

double fd_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept
{
    return fd_from_x(beta * (e_k + (spin_shift * delta) - mu));
}

} // namespace qss::detail
