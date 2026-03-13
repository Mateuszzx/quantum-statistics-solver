/**
 * @file nfl_reworked.cpp
 * @brief Stable implementation of the Non-Fermi-Liquid distribution.
 *
 * This version implements the explicit stationary solution derived from
 *
 *   (1 - n_k)^2 / [ n_{k\sigma} (1 - n_{k\sigma}) ] = exp[ beta (eps_k - mu_sigma) ]
 *
 * with spin competition through n_k = n_{k\uparrow} + n_{k\downarrow} <= 1.
 *
 * For the current interface we keep the convention already used in the code:
 *
 *   mu_sigma = mu - spin_shift * delta,
 *
 * equivalently
 *
 *   eps_{k,sigma} = eps_k + spin_shift * delta.
 */

#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace qss::detail
{

namespace
{
constexpr double LOG_2 = 0.693147180559945309417232121458176568;
constexpr double LOG_4 = 2.0 * LOG_2;
constexpr double LOG_DBL_MIN = -708.396418532264106224411228130256452;

inline double log1p_exp(double x) noexcept
{
    if (x == -std::numeric_limits<double>::infinity())
        return 0.0;

    if (x > 0.0)
        return x + std::log1p(std::exp(-x));

    return std::log1p(std::exp(x));
}

inline double log_cosh(double x) noexcept
{
    const double ax = std::abs(x);

    if (ax < 20.0)
        return std::log(std::cosh(x));

    return ax - LOG_2 + std::log1p(std::exp(-2.0 * ax));
}

inline double log_sinh_abs(double x) noexcept
{
    const double ax = std::abs(x);

    if (ax == 0.0)
        return -std::numeric_limits<double>::infinity();

    if (ax < 20.0)
        return std::log(std::sinh(ax));

    return ax - LOG_2 + std::log1p(-std::exp(-2.0 * ax));
}

inline double nfl_from_stationary_solution(double e_k,
                                           double mu,
                                           double kBT,
                                           double delta,
                                           int spin_shift) noexcept
{
    if (!(kBT > 0.0) || !std::isfinite(kBT))
        return std::numeric_limits<double>::quiet_NaN();

    const double s = static_cast<double>(spin_shift);
    const double a = (mu - e_k) / kBT;
    const double d = -(s * delta) / kBT;
    const double half_d = 0.5 * d;

    const double log_u = LOG_4 + a + 2.0 * log_sinh_abs(half_d);
    const double log_v = LOG_4 + a + 2.0 * log_cosh(half_d);

    const double log_1pu = log1p_exp(log_u);
    const double log_1pv = log1p_exp(log_v);

    const double log_w = -0.5 * (log_1pu + log_1pv);
    const double w = (log_w < LOG_DBL_MIN) ? 0.0 : std::exp(log_w);

    double q = 0.0;
    if (d != 0.0)
    {
        const double log_gu = log_u - log_1pu;
        const double log_gv = log_v - log_1pv;
        const double log_q_abs = 0.5 * (log_gu + log_gv);
        const double q_abs = (log_q_abs < LOG_DBL_MIN) ? 0.0 : std::exp(log_q_abs);
        q = std::copysign(q_abs, d);
    }

    const double n = 0.5 * (1.0 - w + q);
    return std::clamp(n, 0.0, 1.0);
}

} // namespace

double nfl_reworked_distribution(double e_k,
                                 double mu,
                                 double kBT,
                                 double delta,
                                 int spin_shift) noexcept
{
    return nfl_from_stationary_solution(e_k, mu, kBT, delta, spin_shift);
}

} // namespace qss::detail
