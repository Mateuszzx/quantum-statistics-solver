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
/**
 * @brief Core calculation for Non-Fermi-Liquid distribution (algebraic closed form).
 *
 * This implements the MathML solution using auxiliary quantities σ1..σ4.
 *
 * Mapping vs. existing variables:
 *   - MathML uses x, δ (polynomial form)
 *   - Here we set:
 *       x_m = exp(-x)              (same as previous `ex`)
 *       δ_m = exp(x_prim / 2)      so that δ_m^2 = exp(x_prim) (previous `ex_prim`)
 *
 * The MathML yields two algebraic branches (using σ3 and σ2). We select the
 * physically admissible value in [0, 1] when possible.
 */
// inline double nfl_from_x(double x, double x_prim) noexcept
// {
//     if (x > UP_EXP_LIMIT)
//         return 0.0;

//     // Keep your asymptotic guard (unchanged) for extremely negative x.
//     const double ex_prim = std::exp(x_prim);
//     if (x < DOWN_EXP_LIMIT)
//         return std::clamp(0.5 + 1e8 * (1 - ex_prim) / (1 + std::sqrt(ex_prim)), 0.0, 1.0);

//     // --- Map to MathML variables ---
//     const double x_m = std::exp(x);            // x in the MathML-polynomial form
//     const double d_m = std::exp(0.5 * x_prim);  // δ in the MathML-polynomial form

//     // Helper lambdas for powers (avoids repeated pow()).
//     const double d2 = d_m * d_m;
//     const double d3 = d2 * d_m;
//     const double d4 = d2 * d2;
//     const double x2 = x_m * x_m;

//     // Denominator common term: (-δ^2 + x δ + 1)
//     const double denom = (-d2 + x_m * d_m + 1.0);

//     // If denom is ~0, fall back to your previous stable expression.
//     // (This prevents blow-ups on rare pathological parameter combinations.)
//     if (std::abs(denom) < 1e-15)
//     {
//         const double ex = x_m;
//         return std::clamp(ex + 0.5 - std::sqrt(ex_prim * ex * ex + 0.25), 0.0, 1.0);
//     }

//     // σ4 = sqrt( (δ x - 2δ + δ^2 + 1) * (2δ + δ x + δ^2 + 1) )
//     const double s4_arg1 = (d_m * x_m - 2.0 * d_m + d2 + 1.0);
//     const double s4_arg2 = (2.0 * d_m + d_m * x_m + d2 + 1.0);
//     double sigma4_arg = s4_arg1 * s4_arg2;

//     // Numerical safety: tiny negative from rounding -> clamp to 0.
//     if (sigma4_arg < 0.0 && sigma4_arg > -1e-15)
//         sigma4_arg = 0.0;

//     // If genuinely negative, fall back.
//     if (sigma4_arg < 0.0)
//     {
//         const double ex = x_m;
//         return std::clamp(ex + 0.5 - std::sqrt(ex_prim * ex * ex + 0.25), 0.0, 1.0);
//     }

//     const double sigma4 = std::sqrt(sigma4_arg);

//     // σ1 = 2 * (δ^4 + 2 δ^3 x + δ^2 x^2 - 2 δ^2 + 2 δ x + 1)
//     const double sigma1 =
//         2.0 * (d4 + 2.0 * d3 * x_m + d2 * x2 - 2.0 * d2 + 2.0 * d_m * x_m + 1.0);

//     // If σ1 is ~0, fall back.
//     if (std::abs(sigma1) < 1e-15)
//     {
//         const double ex = x_m;
//         return std::clamp(ex + 0.5 - std::sqrt(ex_prim * ex * ex + 0.25), 0.0, 1.0);
//     }

//     // σ2 = δ^2 σ4 - σ4 + 2 δ x + 2 δ^3 x - 2 δ^2 + δ^4 + δ^2 x^2 - δ x σ4 + 1
//     const double sigma2 =
//         d2 * sigma4 - sigma4
//         + 2.0 * d_m * x_m + 2.0 * d3 * x_m
//         - 2.0 * d2 + d4 + d2 * x2
//         - d_m * x_m * sigma4
//         + 1.0;

//     // σ3 = σ4 - δ^2 σ4 + 2 δ x + 2 δ^3 x - 2 δ^2 + δ^4 + δ^2 x^2 + δ x σ4 + 1
//     const double sigma3 =
//         sigma4 - d2 * sigma4
//         + 2.0 * d_m * x_m + 2.0 * d3 * x_m
//         - 2.0 * d2 + d4 + d2 * x2
//         + d_m * x_m * sigma4
//         + 1.0;

//     // Branch formula pattern:
//     // ( (δ^2 * σk)/σ1 - δ^2 - (σk/σ1) + (δ x σk)/σ1 + 1 ) / ( -δ^2 + x δ + 1 )
//     const auto branch_value = [&](double sigmak) -> double {
//         const double sig_over_s1 = sigmak / sigma1;
//         const double num =
//             (d2 * sig_over_s1)
//             - d2
//             - sig_over_s1
//             + (d_m * x_m * sig_over_s1)
//             + 1.0;
//         return num / denom;
//     };

//     const double cand_a = branch_value(sigma3);
//     const double cand_b = branch_value(sigma2);

//     // Select a physically admissible candidate in [0,1] if possible.
//     const bool a_ok = std::isfinite(cand_a) && (cand_a >= 0.0) && (cand_a <= 1.0);
//     const bool b_ok = std::isfinite(cand_b) && (cand_b >= 0.0) && (cand_b <= 1.0);

//     double out = 0.0;
//     if (a_ok && b_ok)
//         out = std::min(cand_a, cand_b);
//     else if (a_ok)
//         out = cand_a;
//     else if (b_ok)
//         out = cand_b;
//     else
//         out = std::min(cand_a, cand_b);

//     return std::clamp(out, 0.0, 1.0);
// }


} // namespace


double nfl_distribution(double e_k, double mu, double kBT, double delta, int spin_shift) noexcept
{
    return nfl_from_x(
        (e_k + (spin_shift * delta) - mu) / kBT,
        (2*(spin_shift * delta)) / kBT
    );
}

} // namespace qss::detail
