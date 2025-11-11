// include/qss/chemical_potential.hpp
#pragma once
#include "qss/distributions.hpp"
#include "qss/solver_params.hpp"
#include "qss/types.hpp"
#include "solver.hpp"
#include <vector>
#include <limits>

namespace qss
{

inline constexpr double QSS_MU_B_EV_PER_T = 5.788381806e-5;

struct ChemicalPotentialResult
{
    double mu = std::numeric_limits<double>::quiet_NaN();
    double error = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool converged = false;
};

// Single beta
ChemicalPotentialResult chemical_potential_B(
    DistType dist,
    double beta,
    double n_assumed,
    const double N,
    const Mat &E,
    double H_a_T,
    double g,
    double muB_eV_per_T,
    const CombinedSolverOptions &solver_opts,
    double a0 = -10.0,
    double b0 = 10.0,
    double x0 = 0.0);

// Vector of betas
std::vector<ChemicalPotentialResult> chemical_potential_B(
    DistType dist,
    const Vec &betas,
    double n_assumed,
    const double N,
    const Mat &E,
    double H_a_T,
    double g,
    double muB_eV_per_T,
    const CombinedSolverOptions &solver_opts,
    double a0 = -10.0,
    double b0 = 10.0,
    double x0 = 0.0);

std::vector<ChemicalPotentialResult> chemical_potential_B_multi_n(
    DistType dist,
    double beta,
    const Vec &n_assumed_values,
    const double N,
    const Mat &E,
    double H_a_T,
    double g,
    double muB_eV_per_T,
    const CombinedSolverOptions &solver_opts,
    double a0 = -10.0,
    double b0 = 10.0,
    double x0 = 0.0);

} // namespace qss
