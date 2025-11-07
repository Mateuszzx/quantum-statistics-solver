#include "../config.hpp" // QSS_OMP_THRESHOLD
#include "impl.hpp"
#include "qss/types.hpp"
#include <cmath>

namespace qss::detail
{

inline constexpr double FD_SPIN = 2.0;

namespace
{
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
#if defined(_OPENMP)
#pragma omp parallel for schedule(static) if (size >= QSS_OMP_THRESHOLD)
#endif
    for (std::size_t i = 0; i < size; ++i)
    {
        out[i] = fd_from_x(beta * (ek[i] - mu));
    }
}

} // namespace qss::detail
