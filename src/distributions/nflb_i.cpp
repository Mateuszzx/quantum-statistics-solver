/**
 * @file nflb_i.cpp
 * @brief Implementation of the Implicit Non-Fermi-Liquid  with Boltzmann(NFLB_i) distribution.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

// ---- Implicit NFLB_i ----
// Solve on f ∈ [0, 1]:
//   p0 f^3 + p1 f^2 + p2 f + p3 = 0
// with p0 = e^x, p1 = 2 e^x + 8, p2 = -16, p3 = 8
// Return 2 * f_root (so result is in [0, 2]).
// Saturation: x >> 0 -> 0, x << 0 -> 2.

namespace
{
/**
 * @brief Finds a root of a cubic polynomial in [0, 1].
 *
 * Solves \f$ ((p_0 f + p_1) f + p_2) f + p_3 = 0 \f$ for \f$ f \in [0, 1] \f$.
 * Uses a robust bisection method.
 *
 * @param p0 Coefficient of \f$ f^3 \f$.
 * @param p1 Coefficient of \f$ f^2 \f$.
 * @param p2 Coefficient of \f$ f^1 \f$.
 * @param p3 Constant term.
 * @param[out] root The found root.
 * @return true if a root was found in [0, 1], false otherwise.
 */
inline bool find_root_on_unit_interval(double p0,
                                       double p1,
                                       double p2,
                                       double p3,
                                       double &root) noexcept
{
    auto poly = [&](double f) noexcept {
        return ((p0 * f + p1) * f + p2) * f + p3;
    };

    double a = 0.0, b = 1.0;
    double fa = poly(a), fb = poly(b);
    if (std::isnan(fa) || std::isnan(fb))
        return false;
    if (fa == 0.0)
    {
        root = a;
        return true;
    }
    if (fb == 0.0)
    {
        root = b;
        return true;
    }
    if (fa * fb > 0.0)
        return false; // no sign change

    // Robust bisection
    for (int it = 0; it < 80; ++it)
    {
        const double m = 0.5 * (a + b);
        const double fm = poly(m);
        if (std::abs(fm) < 1e-12)
        {
            root = m;
            return true;
        }
        if (fa * fm <= 0.0)
        {
            b = m;
            fb = fm;
        }
        else
        {
            a = m;
            fa = fm;
        }
    }
    root = 0.5 * (a + b);
    return true;
}

/**
 * @brief Core calculation for Implicit Non-Fermi-Liquid distribution.
 *
 * Solves the cubic equation:
 * \f[
 * e^x f^3 + (2e^x + 8) f^2 - 16 f + 8 = 0
 * \f]
 * for \f$ f \in [0, 1] \f$, then returns \f$ 2f \f$ (scaled to [0, 2]).
 *
 * @param x Scaled energy difference \f$ \beta(E - \mu) \f$.
 * @return double Distribution value in [0, 2].
 */
inline double nflb_i_from_x(double x) noexcept
{
    if (x > UP_EXP_LIMIT)
        return 0.0; // -> 0
    if (x < DOWN_EXP_LIMIT)
        return 2.0; // -> 2 (spin-like scaling)

    const double ex = std::exp(x);
    const double p0 = ex;
    const double p1 = 2.0 * ex + 8.0;
    const double p2 = -16.0;
    const double p3 = 8.0;

    double froot = 0.0;
    if (find_root_on_unit_interval(p0, p1, p2, p3, froot) && froot >= 0.0 &&
        froot <= 1.0)
    {
        return 2.0 * froot; // scale to [0, 2]
    }
    return 0.0;
}
} // namespace

double nflb_i_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept
{
    return nflb_i_from_x(beta * (e_k + (spin_shift * delta) - mu));
}

} // namespace qss::detail
