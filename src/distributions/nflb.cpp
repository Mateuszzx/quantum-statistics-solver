/**
 * @file nflb.cpp
 * @brief Implementation of the Non-Fermi-Liquid (NFLB) distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

// ---- Non-Fermi-Liquid (NFLB) ----
// For x = beta*(e_k - mu):
//   f(x) = (e^x - sqrt(e^{2x} + 8 e^x) + 8) / (e^x + 8)
// Saturation: x >> 0  -> 0,  x << 0 -> 1

namespace
{
/**
 * @brief Core calculation for Non-Fermi-Liquid distribution.
 * 
 * Computes the distribution function derived from a specific spin-liquid ansatz:
 * \f[
 * f(x) = \frac{e^x - \sqrt{e^{2x} + 8e^x} + 8}{e^x + 8}
 * \f]
 * where \f$ x = \beta(E - \mu) \f$.
 * 
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 1].
 */
inline double nflb_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 1.0;

    const double ex = std::exp(x);
    // Numerically OK up to x≈20 (your clamp); avoids overflow.
    const double num = ex - std::sqrt(ex * ex + 8.0 * ex) + 8.0;
    const double den = ex + 8.0;
    return num / den; // in [0, 1]
}
} // namespace

double nflb_distribution(double e_k, double mu, double beta) noexcept
{
    return nflb_from_x(beta * (e_k - mu));
}

void nflb_distribution_into(const double *ek,
                            std::size_t size,
                            double mu,
                            double beta,
                            double *out)
{
    for (std::size_t i = 0; i < size; ++i)
    {
        out[i] = nflb_from_x(beta * (ek[i] - mu));
    }
}

} // namespace qss::detail
