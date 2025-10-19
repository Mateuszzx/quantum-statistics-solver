#pragma once
#include <optional>
#include <tuple>
#include "types.hpp"

namespace qss {

enum class DistType { FD, FLB, NFL, NFL_i };

// --- Core scalar/vector distributions (public API) ---
double fd_distribution(double e_k, double mu, double beta);
Vec    fd_distribution(const Vec& e_k, double mu, double beta);

double flb_distribution(double e_k, double mu, double beta);
Vec    flb_distribution(const Vec& e_k, double mu, double beta);

double nfl_distribution(double e_k, double mu, double beta);
Vec    nfl_distribution(const Vec& e_k, double mu, double beta);

double nfl_i_distribution(double e_k, double mu, double beta);
Vec    nfl_i_distribution(const Vec& e_k, double mu, double beta);

// --- Dispatcher helpers (public, optional) ---
std::optional<double> distribution(double e_k, double mu, double beta, DistType dist);
std::optional<Vec>    distribution(const Vec& e_k, double mu, double beta, DistType dist);

// --- 2D symmetry helper (public) ---
Mat distribution_symmetry_2d(const Mat& e_k, double mu, double beta, DistType dist);

// --- Lattice preparation (public) ---
std::tuple<Mat, Mat, Mat> prepare_2D_lattice_data(const Vec& kx, const Vec& ky,
                                                  double t, double t_prim);

} // namespace qss
