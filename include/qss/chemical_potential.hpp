// include/qss/chemical_potential.hpp
#pragma once
#include "qss/types.hpp"
#include "qss/distributions.hpp"
#include "qss/solver_params.hpp"

namespace qss {

// Bohr magneton in eV/T (default units consistent with E grid in eV)
inline constexpr double QSS_MU_B_EV_PER_T = 5.788381806e-5;

// 2D only; H_a in Tesla; Zeeman computed inside: Δ = 0.5 * g * μ_B * H_a
double chemical_potential_B(
    DistType dist,
    double beta,
    double n_assumed,
    const Vec& kx, const Vec& ky,
    const Mat& E,
    double H_a_T,                         // magnetic field (Tesla)
    double g = 2.0,                       // g-factor
    double muB_eV_per_T = QSS_MU_B_EV_PER_T,
    const SolverParams& params = {});

// Vectorized over betas (same H_a for all)
Vec chemical_potential_B(
    DistType dist,
    const Vec& betas,
    double n_assumed,
    const Vec& kx, const Vec& ky,
    const Mat& E,
    double H_a_T,
    double g = 2.0,
    double muB_eV_per_T = QSS_MU_B_EV_PER_T,
    const SolverParams& params = {});

} // namespace qss
