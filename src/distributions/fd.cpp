/**
 * @file fd.cpp
 * @brief Implementation of the Fermi-Dirac distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

/// Spin degeneracy factor for Fermi-Dirac statistics.
inline constexpr double FD_SPIN = 2.0;

namespace
{
/**
 * @brief Core calculation for Fermi-Dirac distribution.
 * 
 * Computes \f$ f(x) = \frac{2}{e^x + 1} \f$ where \f$ x = \beta(E - \mu) \f$.
 * Handles numerical stability for large positive and negative x.
 * 
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 2].
 */
inline double fd_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return FD_SPIN;
    if (x >= 0.0)
        return FD_SPIN / (1.0 + std::exp(x));
    const double ex = std::exp(-x);
    return FD_SPIN * ex / (1.0 + ex);
}
} // namespace

double fd_distribution(double e_k, double mu, double beta) noexcept
{
    return fd_from_x(beta * (e_k - mu));
}

void fd_distribution_into(const double *ek,
                          std::size_t size,
                          double mu,
                          double beta,
                          double *out)
{
    for (std::size_t i = 0; i < size; ++i)
    {
        out[i] = fd_from_x(beta * (ek[i] - mu));
    }
}

} // namespace qss::detail
