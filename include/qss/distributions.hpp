#pragma once
#include "types.hpp"
#include <optional>
#include <tuple>

namespace qss
{

enum class DistType
{
    FD,
    FLB,
    NFLB,
    NFLB_i
};

std::optional<double> distribution(double e_k,
                                   double mu,
                                   double beta,
                                   DistType dist) noexcept;

std::optional<Vec> distribution(const Vec &e_k,
                                double mu,
                                double beta,
                                DistType dist);
} // namespace qss
