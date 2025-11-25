/**
 * @file flb.cpp
 * @brief Implementation of the Fermi-Liquid-Like (Spin Liquid) distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

/// Spin degeneracy factor.
inline constexpr double FD_SPIN = 2.0;

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
        return FD_SPIN;
    if (x >= 0.0)
        return FD_SPIN * 2.0 / (2.0 + std::exp(x));
    const double ex = 2.0 * std::exp(-x);
    return FD_SPIN * ex / (1 + ex);
}
} // namespace

double flb_distribution(double e_k, double mu, double beta) noexcept
{
    return flb_from_x(beta * (e_k - mu));
}

void flb_distribution_into(const double *ek,
                           std::size_t size,
                           double mu,
                           double beta,
                           double *out)
{
    for (std::size_t i = 0; i < size; ++i)
    {
        out[i] = flb_from_x(beta * (ek[i] - mu));
    }
}

} // namespace qss::detail
