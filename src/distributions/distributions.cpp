#include "qss/distributions.hpp"
#include "impl.hpp"


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

} // namespace qss
