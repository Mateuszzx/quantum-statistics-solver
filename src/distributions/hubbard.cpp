/**
 * @file hubbard.cpp
 * @brief Implementation of Hubbard distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{


namespace
{
/**
 * @brief Core calculation for Hubbard distribution.
 *
 * Computes \f$ f(x) = \frac{e^x}{e^x + e^x'prim + e^{-U}} \f$ where \f$ x = \beta(E - \mu) \f$.
 * Handles numerical stability for large positive and negative x.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 0.5].
 */
inline double hubbard_from_x(double x, double x_prim, double U, bool n_kd = false) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    if (x < DOWN_EXP_LIMIT)
        return 1. / (1. + std::exp(x-x_prim));

    const double ex = std::exp(-x);
    const double ex_prim = std::exp(-x_prim);
    const double ex_U = ex* ex_prim * std::exp(-U);

    if(!n_kd)
        return (ex + ex_U)/ (1.0 + ex + ex_prim + ex_U);
    else
        return ex_U/ (1.0 + ex + ex_prim + ex_U);
}
} // namespace

double hubbard_distribution(double e_k, double mu, double beta, double delta, int spin_shift, double U, bool n_kd) noexcept
{
    return hubbard_from_x(beta * (e_k + (spin_shift * delta) - mu), beta * (e_k - (spin_shift * delta) - mu), U, n_kd);
}

} // namespace qss::detail
