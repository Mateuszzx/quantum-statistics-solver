#pragma once
#include <optional>
#include <tuple>
#include "types.hpp"

namespace qss {

enum class DistType { FD, FLB, NFLB, NFLB_i };

std::optional<double> distribution(double e_k, double mu, double beta, DistType dist) noexcept;
std::optional<Vec>    distribution(const Vec& e_k, double mu, double beta, DistType dist);

// Keep your lattice helper public if you want:
std::tuple<Mat, Mat, Mat> prepare_2D_lattice_data(const Vec& kx, const Vec& ky,
                                                  double t, double t_prim);

} // namespace qss
