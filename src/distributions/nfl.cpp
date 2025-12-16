/**
 * @file nfl.cpp
 * @brief Implementation of the Non-Fermi-Liquid distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <algorithm>
#include <cmath>

namespace qss::detail
{

namespace
{
/**
 * @brief Core calculation for Non-Fermi-Liquid distribution.
 *
 * Computes \f$ f(x) = \frac{1/2} + e^-x - 0.5* \sqrt{1 + \frac{4}{e^{-x} e{-x_prim}}}\f$ where \f$ x = \beta(E - \mu) \f$.
 * This distribution arises in certain spin liquid models.
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 0.5].
 */
inline double nfl_from_x(double x, double x_prim) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0;
    const double ex_prim = std::exp(x_prim);
    if (x < DOWN_EXP_LIMIT)
        return std::clamp(0.5 +  1e8 * (1-ex_prim)/(1+std::sqrt(ex_prim)), 0.0, 1.0);

    const double ex = std::exp(-x);


    return std::clamp(ex +  0.5 - std::sqrt(ex_prim * ex* ex + 0.25), 0.0, 1.0);
}
} // namespace


double nfl_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept
{
    return nfl_from_x(
        beta * (e_k + (spin_shift * delta) - mu),
        beta * (2*(spin_shift * delta))
    );
}

} // namespace qss::detail
