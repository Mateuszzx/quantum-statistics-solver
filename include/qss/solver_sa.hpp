#pragma once
#include "qss/solver_types.hpp"
#include <string>

namespace qss
{


/**
 * @brief Cooling schedule for Simulated Annealing.
 */
enum class CoolingSchedule
{
    Exponential, ///< T = T0 * alpha^k
    Linear,      ///< T = T0 * (1 - k/max_iters)
    Logarithmic, ///< T = T0 / log(2 + k)
    Sawtooth     ///< Periodic linear decay
};

// Derive from SolverOptions (fix syntax)
/**
 * @brief Configuration options for Simulated Annealing optimization.
 */
struct SimulatedAnnealingOptions : public SolverOptions
{
    SimulatedAnnealingOptions()
    {
        // Set default iterations for SA (higher than generic root-finding default)
        max_iters = 50'000;
    }

    double T0 = 1.0;        ///< Initial temperature.
    double T_min = 1e-6;    ///< Minimum temperature to stop.
    double alpha = 0.95;    ///< Decay factor for Exponential schedule.
    CoolingSchedule sched = CoolingSchedule::Exponential; ///< Cooling schedule type.

    std::size_t max_iters = 50'000;  ///< Maximum total iterations.
    std::size_t min_iters = 1'000;   ///< Minimum iterations before convergence checks.
    std::size_t stall_iters = 5'000; ///< Stop if no improvement after this many iterations.

    double step_sigma = 0.1;      ///< Base standard deviation for proposal step.
    double step_decay = 1.0;      ///< Decay factor for step size per iteration (<= 1.0).
    bool scale_sigma_by_T = true; ///< If true, scale step size by sqrt(T/T0).

    std::vector<std::pair<double, double>>
        bounds; // size=D or empty; +/-inf allowed  ///< Bounds for each dimension [min, max].
    bool reflect_bounds = true; ///< If true, reflect off bounds; otherwise clamp.

    unsigned long long seed = 0xC0FFEEULL; ///< Random seed.

    std::string log_file = ""; ///< Path to log file. If empty, logging is disabled.
};

// Declare SA entry point used by OptimizationSolver
/**
 * @brief Internal implementation of Simulated Annealing minimization.
 *
 * @param g The objective function.
 * @param x0 Initial guess.
 * @param o SA options.
 * @return MinResult Optimization result.
 */
MinResult sa_minimize_impl(const ObjectiveND &g,
                           const std::vector<double> &x0,
                           const SimulatedAnnealingOptions &o);

} // namespace qss
