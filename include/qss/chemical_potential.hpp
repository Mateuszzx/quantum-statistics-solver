// include/qss/chemical_potential.hpp
#pragma once
#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include "solver.hpp"
#include <vector>
#include <limits>

namespace qss
{

inline constexpr double QSS_MU_B_EV_PER_T = 5.788381806e-5;

/**
 * @brief Result structure for chemical potential calculations.
 */
struct ChemicalPotentialResult
{
    double mu = std::numeric_limits<double>::quiet_NaN(); ///< The calculated chemical potential.
    double error = std::numeric_limits<double>::quiet_NaN(); ///< Absolute error |n_calc - n_target|.
    int iterations = 0; ///< Number of iterations performed by the solver.
    bool converged = false; ///< True if the solver converged successfully.
};

/**
 * @brief Calculate chemical potential for a single set of parameters.
 *
 * Solves for mu such that the total particle density equals n_assumed.
 *
 * @param dist Distribution type (FD, FLB, NFLB, etc.).
 * @param beta Inverse temperature (1/k_B T).
 * @param n_assumed Target particle density (filling factor).
 * @param N Linear size of the lattice (used for validation).
 * @param E Energy dispersion matrix (2D).
 * @param H_a_T Applied magnetic field in Tesla.
 * @param g Landé g-factor.
 * @param muB_eV_per_T Bohr magneton in eV/T.
 * @param solver_opts Configuration options for the numerical solver.
 * @param a0 Lower bound for root finding (default -10.0).
 * @param b0 Upper bound for root finding (default 10.0).
 * @param x0 Initial guess for minimization (default 0.0).
 * @return ChemicalPotentialResult Result containing mu and status.
 */
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


/**
 * @brief Calculate chemical potential for a vector of inverse temperatures (betas).
 *
 * @param dist Distribution type.
 * @param betas Vector of inverse temperatures.
 * @param n_assumed Target particle density.
 * @param N Linear size of the lattice.
 * @param E Energy dispersion matrix.
 * @param H_a_T Applied magnetic field in Tesla.
 * @param g Landé g-factor.
 * @param muB_eV_per_T Bohr magneton in eV/T.
 * @param solver_opts Solver options.
 * @param a0 Lower bound.
 * @param b0 Upper bound.
 * @param x0 Initial guess.
 * @return std::vector<ChemicalPotentialResult> Vector of results corresponding to betas.
 */
std::vector<ChemicalPotentialResult> chemical_potential_B_multi_beta(
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

/**
 * @brief Calculate chemical potential for a vector of target densities (n).
 *
 * @param dist Distribution type.
 * @param beta Inverse temperature.
 * @param n_assumed_values Vector of target particle densities.
 * @param N Linear size of the lattice.
 * @param E Energy dispersion matrix.
 * @param H_a_T Applied magnetic field in Tesla.
 * @param g Landé g-factor.
 * @param muB_eV_per_T Bohr magneton in eV/T.
 * @param solver_opts Solver options.
 * @param a0 Lower bound.
 * @param b0 Upper bound.
 * @param x0 Initial guess.
 * @return std::vector<ChemicalPotentialResult> Vector of results corresponding to n_assumed_values.
 */
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
