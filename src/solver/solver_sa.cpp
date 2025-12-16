#include "qss/solver_sa.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <stdexcept>

namespace qss
{

/**
 * @brief Reflects a value x back into the interval [lo, hi].
 *
 * This function implements a reflection logic. If x is outside [lo, hi],
 * it is reflected back into the domain. This preserves the detailed balance condition better
 * than simple clamping (which accumulates probability mass at the boundaries).
 *
 * @param x The value to reflect.
 * @param lo The lower bound.
 * @param hi The upper bound.
 * @return The reflected value within [lo, hi].
 */
static inline double reflect_box(double x, double lo, double hi)
{
    if (!(std::isfinite(lo) && std::isfinite(hi)) || lo >= hi)
        return x;
    const double L = hi - lo;
    double y = std::fmod(x - lo, 2 * L);
    if (y < 0)
        y += 2 * L;
    return (y <= L) ? (lo + y) : (hi - (y - L));
}

/**
 * @brief Calculates the temperature for the k-th iteration based on the cooling schedule.
 *
 * @param o The simulated annealing options containing schedule parameters.
 * @param k The current iteration index (0-based).
 * @return The temperature T for the current step.
 */
static inline double next_T(const SimulatedAnnealingOptions &o, std::size_t k)
{
    switch (o.sched)
    {
    case CoolingSchedule::Exponential:
        // Classic geometric cooling: T(k) = T0 * alpha^k
        return o.T0 * std::pow(o.alpha, (double)k);
    case CoolingSchedule::Linear:
    {
        // Linear decrease: T(k) = T0 * (1 - k/max_iters)
        if (o.max_iters == 0)
            return 0.0;
        double frac = 1.0 - (double)k / (o.alpha * (double)o.max_iters);
        return std::max(0.0, o.T0 * frac);
    }
    case CoolingSchedule::Logarithmic:
        // Geman & Geman schedule: T(k) = T0 / log(k + 2)
        // Theoretically guarantees convergence to global optimum if T0 is high enough,
        // but is very slow in practice.
        return o.T0 / std::log(2.0 + (double)k/o.alpha);
    case CoolingSchedule::Sawtooth:
    {
        // Periodic reheating (Restart strategy).
        // Useful to escape deep local minima by periodically raising T.
        if (o.max_iters == 0)
            return 0.0;
        // Fixed 5 cycles for now
        const double cycles = 5.0;
        double cycle_len = (double)o.max_iters / cycles;
        double rel_k = std::fmod((double)k, cycle_len);
        return o.T0 * (1.0 - rel_k / (o.alpha*cycle_len));
    }
    }
    return 0.0;
}

/**
 * @brief Core implementation of the Simulated Annealing minimization algorithm.
 *
 * This function attempts to find the global minimum of an N-dimensional objective function
 * using a probabilistic technique that approximates the global optimum.
 *
 * Algorithm Overview:
 * 1. Initialization: Starts with an initial solution x0 and an initial temperature T0.
 * 2. Iteration:
 *    - Temperature Update: The temperature T decreases according to a cooling schedule.
 *      As T approaches zero, the algorithm becomes more greedy.
 *    - Step Size Adaptation: The search radius (sigma) decays over time to refine the search.
 *    - Candidate Generation: A new candidate is proposed by perturbing the current solution
 *      with random noise. Boundary conditions are handled (clamping or reflection).
 *    - Evaluation: The objective function is evaluated at the candidate point.
 *    - Acceptance Criteria (Metropolis):
 *      - If the candidate is better (lower energy), it is always accepted.
 *      - If the candidate is worse, it is accepted with probability P = exp(-(E_new - E_old) / T).
 *        This allows escaping local minima.
 * 3. Termination: Stops on max iterations, min temperature, or stall iterations.
 *
 * @param g The objective function to minimize.
 * @param x0 The starting point.
 * @param o Configuration options.
 * @return MinResult Structure containing the best found solution and statistics.
 */
MinResult sa_minimize_impl(const ObjectiveND &g,
                           const std::vector<double> &x0,
                           const SimulatedAnnealingOptions &o)
{
    MinResult R{};
    if (x0.empty())
    {
        R.status = Status::DomainError;
        return R;
    }

    std::mt19937_64 rng((std::mt19937_64::result_type)o.seed);
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::normal_distribution<double> norm(0.0, 1.0);

    // Initialize best solution with the starting point
    R.x_best = x0;
    double f_best = g(R.x_best);
    if (!std::isfinite(f_best))
    {
        R.status = Status::DomainError;
        return R;
    }
    R.f_best = f_best;

    // Current working solution
    std::vector<double> x_cur = R.x_best;
    double f_cur = f_best;

    std::vector<double> x_cand(x0.size());
    std::size_t no_improve = 0;
    double T = o.T0;
    const double base_sigma = std::max(1e-16, o.step_sigma);
    double sigma = base_sigma;

    // Setup logging if requested
    std::ofstream log_stream;
    if (!o.log_file.empty())
    {
        log_stream.open(o.log_file);
        if (log_stream.is_open())
        {
            log_stream << "iteration,temperature,energy";
            for (size_t i = 0; i < x0.size(); ++i)
            {
                log_stream << ",x" << i;
            }
            log_stream << "\n";
        }
    }

    // Main Optimization Loop
    // Cast max_iters to size_t to avoid signed/unsigned comparison warnings
    const std::size_t max_iters = static_cast<std::size_t>(o.max_iters);
    for (std::size_t k = 0; k < max_iters; ++k)
    {
        R.iterations = k + 1;

        // --- Step 1: Update Temperature ---
        // The temperature controls the probability of accepting worse solutions.
        // High T -> High probability of accepting uphill moves (exploration).
        // Low T -> Low probability (exploitation).
        T = next_T(o, k);

        // Log current state
        if (log_stream.is_open())
        {
            log_stream << (k + 1) << "," << T << "," << f_cur;
            for (double val : x_cur)
            {
                log_stream << "," << val;
            }
            log_stream << "\n";
        }

        // --- Step 2: Update Step Size (Sigma) ---
        // The step size determines how far we look for a neighbor.
        // It decays over time to allow fine-tuning of the solution in later stages.
        sigma = base_sigma;
        if (o.scale_sigma_by_T && o.T0 > 0.0)
        {
            sigma *= std::sqrt(std::max(0.0, T / o.T0));
        }

        const double decay_k = std::pow(std::max(0.0, o.step_decay), static_cast<double>(k));
        sigma *= decay_k;
        sigma = std::max(1e-8, sigma);

        if (k + 1 >= o.min_iters)
        {
            // For Sawtooth, T drops to 0 periodically. We shouldn't stop early on T_min.
            if (o.sched != CoolingSchedule::Sawtooth && T < o.T_min)
                break;
            if (no_improve >= o.stall_iters)
                break;
        }

        // --- Step 3: Propose Candidate Solution ---
        // Perturb the current solution x_cur by adding random noise (Gaussian).
        // Also enforce boundary conditions (reflect or clamp).
        for (std::size_t d = 0; d < x_cur.size(); ++d)
        {
            double step = sigma * norm(rng);
            double xc = x_cur[d] + step;
            if (!o.bounds.empty())
            {
                auto [lo, hi] = o.bounds[d];
                if (o.reflect_bounds && std::isfinite(lo) && std::isfinite(hi))
                    xc = reflect_box(xc, lo, hi);
                else
                {
                    if (std::isfinite(lo) && xc < lo)
                        xc = lo;
                    if (std::isfinite(hi) && xc > hi)
                        xc = hi;
                }
            }
            x_cand[d] = xc;
        }

        double f_cand = g(x_cand);
        if (!std::isfinite(f_cand))
        {
            ++R.rejects;
            continue;
        }

        // --- Step 4: Metropolis Acceptance Criterion ---
        // Calculate energy difference dE = f_cand - f_cur.
        // - If dE <= 0 (improvement), we always accept.
        // - If dE > 0 (worse), we accept with probability exp(-dE / T).
        double dE = f_cand - f_cur;
        bool accept =
            (dE <= 0.0) || (uni(rng) < std::exp(-dE / std::max(T, 1e-8)));

        if (accept)
        {
            std::swap(x_cur, x_cand);
            f_cur = f_cand;
            ++R.accepts;
            if (f_cur < f_best)
            {
                // Check for significant improvement using SolverOptions tolerances.
                // If the improvement is smaller than tolerance, we don't reset the stall counter.
                // This prevents the solver from running forever on tiny improvements (e.g. numerical noise).
                double improvement = f_best - f_cur;
                double tol = o.abs_tol + o.rel_tol * std::abs(f_best);

                f_best = f_cur;
                R.x_best = x_cur;
                R.f_best = f_best;

                if (improvement > tol)
                    no_improve = 0;
                else
                    ++no_improve;
            }
            else
            {
                ++no_improve;
            }
        }
        else
        {
            ++R.rejects;
            ++no_improve;
        }
    }

    R.T_final = T;
    R.status = Status::Converged; // or differentiate by stop reason if you want
    return R;
}

} // namespace qss
