/**
 * @file distributions.cpp
 * @brief Dispatch logic for statistical distributions.
 * 
 * Implements the main `distribution` function overloads that dispatch to specific
 * distribution implementations (FD, FLB, NFLB, NFLB_i) based on the `DistType` enum.
 */

#include "qss/distributions.hpp"


namespace qss
{

std::optional<double> distribution(double e_k,
                                double mu,
                                double beta,
                                DistType dist) noexcept
{
    using namespace qss::detail;
    switch (dist)
    {
    case DistType::FD:
        return fd_distribution(e_k, mu, beta);
    case DistType::FLB:
        return flb_distribution(e_k, mu, beta);
    case DistType::NFLB:
        return nflb_distribution(e_k, mu, beta);
    case DistType::NFLB_i:
        return nflb_i_distribution(e_k, mu, beta);
    default:
        return std::nullopt;
    }
}

std::optional<Vec> distribution(const Vec &e_k,
                                double mu,
                                double beta,
                                DistType dist)
{
    using namespace qss::detail;
    const std::size_t n = e_k.size();
    Vec out(n);
    switch (dist)
    {
    case DistType::FD:
        fd_distribution_into(e_k.data(), n, mu, beta, out.data());
        break;
    case DistType::FLB:
        flb_distribution_into(e_k.data(), n, mu, beta, out.data());
        break;
    case DistType::NFLB:
        nflb_distribution_into(e_k.data(), n, mu, beta, out.data());
        break;
    case DistType::NFLB_i:
        nflb_i_distribution_into(e_k.data(), n, mu, beta, out.data());
        break;
    default:
        return std::nullopt;
    }
    return out;
}

std::optional<Mat> distribution(const Mat &e_k,
                                double mu,
                                double beta,
                                DistType dist)
{
    const std::size_t rows = e_k.size();
    Mat out(rows);
    for (std::size_t i = 0; i < rows; ++i)
    {
        auto row_res = distribution(e_k[i], mu, beta, dist);
        if (!row_res)
            return std::nullopt;
        out[i] = std::move(*row_res);
    }
    return out;
}

} // namespace qss
