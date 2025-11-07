// include/qss/chemical_potential.hpp
#pragma once
#include "qss/distributions.hpp"
#include "qss/solver_params.hpp"
#include "qss/types.hpp"
#include "solver.hpp"

namespace qss
{

// Bohr magneton in eV/T (default units consistent with E grid in eV)
inline constexpr double QSS_MU_B_EV_PER_T = 5.788381806e-5;


double chemical_potential_B(DistType dist,
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

Vec chemical_potential_B(DistType dist,
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

} // namespace qss
