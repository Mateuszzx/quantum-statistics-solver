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
 * Computes occupation number \f$ n_{\sigma} = \frac{e^{-\beta(E_\sigma - \mu)} + e^{-\beta(E_\uparrow + E_\downarrow + U - 2\mu)}}{Z} \f$
 * or double occupation \f$ d = \frac{e^{-\beta(E_\uparrow + E_\downarrow + U - 2\mu)}}{Z} \f$.
 *
 * @param x Scaled energy difference \f$ \beta(E_\sigma - \mu) \f$.
 * @param x_prim Scaled energy difference \f$ \beta(E_{-\sigma} - \mu) \f$.
 * @param u Scaled interaction parameter \f$ \beta U \f$.
 * @param n_kd If true, returns double occupation probability.
 * @return double Distribution value in [0, 1].
 */
inline double hubbard_from_x(double x, double x_prim, double u, bool n_kd = false) noexcept
{
    // Exponents: 0, -x, -x_prim, -x - x_prim - u
    // We shift all exponents by subtracting max_val to avoid overflow.
    // max_val = max(0, -x, -x_prim, -x - x_prim - u)

    double exp1 = -x;
    double exp2 = -x_prim;
    double exp3 = -x - x_prim - u;

    double max_val = 0.0; // Corresponds to the '1' term (exp(0))
    if (exp1 > max_val) max_val = exp1;
    if (exp2 > max_val) max_val = exp2;
    if (exp3 > max_val) max_val = exp3;

    // Calculate terms relative to max_val
    double term0 = std::exp(0.0 - max_val);      // 1
    double term1 = std::exp(exp1 - max_val);     // e^{-x}
    double term2 = std::exp(exp2 - max_val);     // e^{-x_prim}
    double term3 = std::exp(exp3 - max_val);     // e^{-x - x_prim - u}

    double Z = term0 + term1 + term2 + term3;

    if (n_kd)
    {
        return term3 / Z;
    }
    else
    {
        return (term1 + term3) / Z;
    }
}
} // namespace

double hubbard_distribution(double e_k, double mu, double kBT, double delta, int spin_shift, double U, bool n_kd) noexcept
{
    // Pass U/kBT as the scaled interaction parameter
    return hubbard_from_x((e_k + (spin_shift * delta) - mu) / kBT, (e_k - (spin_shift * delta) - mu) / kBT, U / (2*kBT), n_kd);
}

} // namespace qss::detail
