/**
 * @file nflb.cpp
 * @brief Implementation of the Non-Fermi-Liquid with Boltzmann (NFLB) distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

// ---- Non-Fermi-Liquid (NFLB) ----
// For x = beta*(e_k - mu):
//   f(x) = (e^x - sqrt(e^{2x} + 8 e^x) + 8) / (2e^x + 166)
// Saturation: x >> 0  -> 0,  x << 0 -> 1/2

namespace
{
/**
 * @brief Core calculation for Non-Fermi-Liquid distribution.
 *
 * Computes the distribution function derived from a specific spin-liquid ansatz:
 * \f[
 * f(x) = \frac{e^x - \sqrt{e^{2x} + 8e^x} + 8}{2e^x + 16}
 * \f]
 * where \f$ x = \beta(E - \mu) \f$.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 0.5].
 */
inline double nflb_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 0.5;

    const double ex = std::exp(x);
    const double num = ex - std::sqrt(ex * ex + 8.0 * ex) + 8.0;
    const double den = ex + 8.0;
    return num / (2*den); // in [0, 0.5]
}
} // namespace

double nflb_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept
{
    return nflb_from_x(beta * (e_k + (spin_shift * delta) - mu));
}

} // namespace qss::detail

