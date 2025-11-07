#pragma once
#include "solver_types.hpp"

namespace qss
{


enum class CoolingSchedule
{
    Exponential,
    Linear,
    Logarithmic
};

// Derive from SolverOptions (fix syntax)
struct SimulatedAnnealingOptions : public SolverOptions
{
    double T0 = 1.0;
    double T_min = 1e-6;
    double alpha = 0.95;
    CoolingSchedule sched = CoolingSchedule::Exponential;

    std::size_t max_iters = 50'000;
    std::size_t min_iters = 1'000;
    std::size_t stall_iters = 5'000;

    double step_sigma = 0.1;
    double step_decay = 1.0; // <=1 to shrink over time
    bool scale_sigma_by_T = true;

    std::vector<std::pair<double, double>>
        bounds; // size=D or empty; +/-inf allowed
    bool reflect_bounds = true;

    unsigned long long seed = 0xC0FFEEULL;
};

// Declare SA entry point used by OptimizationSolver
MinResult sa_minimize_impl(const ObjectiveND &g,
                           const std::vector<double> &x0,
                           const SimulatedAnnealingOptions &o);

} // namespace qss
