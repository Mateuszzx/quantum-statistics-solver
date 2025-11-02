// include/qss/chemical_potential.hpp
#pragma once
#include "qss/types.hpp"
#include "qss/distributions.hpp"
#include "qss/solver_params.hpp"

namespace qss {

// 2D only; delta = Zeeman split (same units as E)
double chemical_potential(
    DistType dist,
    double beta,
    double n_assumed,
    const Vec& kx, const Vec& ky,
    const Mat& E,
    double delta = 0.0,
    const SolverParams& params = {});

// Vectorized over betas
Vec chemical_potential(
    DistType dist,
    const Vec& betas,
    double n_assumed,
    const Vec& kx, const Vec& ky,
    const Mat& E,
    double delta = 0.0,
    const SolverParams& params = {});

} // namespace qss
