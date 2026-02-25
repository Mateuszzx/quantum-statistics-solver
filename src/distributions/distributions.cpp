/**
 * @file distributions.cpp
 * @brief Dispatch logic for statistical distributions.
 *
 * Implements the main `distribution` function overloads that dispatch to specific
 * distribution implementations (FD, FLB, NFL, NFLB, SSL) based on the `DistType` enum.
 */

#include "qss/distributions.hpp"


namespace qss
{

std::optional<double> distribution(double e_k,
                                double mu,
                                double beta,
                                double delta,
                                int spin_shift,
                                double U,
                                DistType dist,
                                bool n_kd) noexcept
{
    using namespace qss::detail;
    switch (dist)
    {
    case DistType::FD:
        return fd_distribution(e_k, mu, beta, delta, spin_shift);
    case DistType::FLB:
        return flb_distribution(e_k, mu, beta, delta, spin_shift);
    case DistType::NFL:
        return nfl_distribution(e_k, mu, beta, delta, spin_shift);
    case DistType::NFLB:
        return nflb_distribution(e_k, mu, beta, delta, spin_shift);
    case DistType::SSL:
        return ssl_distribution(e_k, mu, beta, delta, spin_shift);
    case DistType::HUBBARD:
        return hubbard_distribution(e_k, mu, beta, delta, spin_shift, U, n_kd);
    default:
        return std::nullopt;
    }
}

std::optional<Vec> distribution(const Vec &e_k,
                                double mu,
                                double beta,
                                double delta,
                                int spin_shift,
                                double U,
                                DistType dist,
                                bool n_kd)
{

    using namespace qss::detail;
    const std::size_t n = e_k.size();
    Vec out(n);

    for (std::size_t i = 0; i < n; ++i)
    {
        auto val = distribution(e_k[i], mu, beta, delta, spin_shift, U, dist, n_kd);
        if (!val)
            return std::nullopt;
        out[i] = *val;
    }

    return out;
}

std::optional<Mat> distribution(const Mat &e_k,
                                double mu,
                                double beta,
                                double delta,
                                int spin_shift,
                                double U,
                                DistType dist,
                                bool n_kd)
{
    const std::size_t rows = e_k.size();
    Mat out(rows);
    for (std::size_t i = 0; i < rows; ++i)
    {
        auto row_res = distribution(e_k[i], mu, beta, delta, spin_shift, U, dist, n_kd);
        if (!row_res)
            return std::nullopt;
        out[i] = std::move(*row_res);
    }
    return out;
}

} // namespace qss
