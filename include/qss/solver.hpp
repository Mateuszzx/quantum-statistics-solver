#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "qss/bracket.hpp"
#include "qss/solver_sa.hpp"

namespace qss
{

// ---------- enums ----------
/**
 * @brief Available solver algorithms.
 */
enum class SolverType
{
    Bisection,
    Secant,
    Brent,
    SimulatedAnnealing
};

// NOTE: in existing code base SolverOptions is used for root-finding options.
// If you already have a specific RootOptions type elsewhere, you can alias it here.
// using RootSolverOptions = RootOptions;
using RootSolverOptions =
    SolverOptions; // keep backward-compat: your existing "SolverOptions"

/**
 * @brief Combined configuration for all solver types.
 */
struct CombinedSolverOptions
{
    SolverType solver = SolverType::Brent;              ///< The solver algorithm to use.
    BracketerType bracketer_type = BracketerType::Expansion; ///< Strategy for bracketing roots.
    BracketerOptions bracketer_opts{};                  ///< Options for the bracketer.
    RootSolverOptions root{};                           ///< Options for root-finding algorithms.
    SimulatedAnnealingOptions sa{};                     ///< Options for Simulated Annealing.
};

// ---------- OptimizationSolver (public interface) ----------
/**
 * @brief Main solver class providing a unified interface for root finding and minimization.
 */
class OptimizationSolver
{
public:
    /**
     * @brief Construct a new Optimization Solver.
     *
     * @param opts Solver configuration options.
     */
    explicit OptimizationSolver(const CombinedSolverOptions &opts = {})
        : opts_(opts), bracketer_(opts_.bracketer_type, opts_.bracketer_opts)
    {
    }

    // configuration
    /**
     * @brief Update the full set of options.
     *
     * @param o New options.
     */
    void set_options(const CombinedSolverOptions &o)
    {
        opts_ = o;
        bracketer_ = Bracketer(opts_.bracketer_type, opts_.bracketer_opts);
    }

    /**
     * @brief Set the solver algorithm type.
     *
     * @param t Solver type.
     */
    void set_solver_type(SolverType t) noexcept
    {
        opts_.solver = t;
    }

    /**
     * @brief Configure the bracketer.
     *
     * @param t Bracketer type.
     * @param o Bracketer options.
     */
    void set_bracketer(BracketerType t, const BracketerOptions &o = {}) noexcept
    {
        opts_.bracketer_type = t;
        opts_.bracketer_opts = o;
        bracketer_ = Bracketer(t, o);
    }

    /**
     * @brief Set options specific to root finding.
     *
     * @param o Root solver options.
     */
    void set_root_options(const RootSolverOptions &o) noexcept
    {
        opts_.root = o;
    }

    /**
     * @brief Set options specific to Simulated Annealing.
     *
     * @param o SA options.
     */
    void set_sa_options(const SimulatedAnnealingOptions &o) noexcept
    {
        opts_.sa = o;
    }

    /**
     * @brief Get the current options.
     *
     * @return const CombinedSolverOptions&
     */
    const CombinedSolverOptions &options() const noexcept
    {
        return opts_;
    }

    /**
     * @brief Get the internal bracketer instance.
     *
     * @return const Bracketer&
     */
    const Bracketer &bracketer() const noexcept
    {
        return bracketer_;
    }

    // ROOT-FINDING (bisection/secant/brent)
    /**
     * @brief Solve for a root f(x) = 0 in the interval [a0, b0].
     *
     * @param f The function to solve.
     * @param a0 Initial lower bound.
     * @param b0 Initial upper bound.
     * @return RootResult Result containing the root and status.
     */
    RootResult solve_root(const Function1D &f, double a0, double b0) const;

    // Unified scalar solve:
    // - if solver == SimulatedAnnealing -> minimize |f(x)| starting at x0 and return x*
    // - else -> root-find f(x)=0 on [a0,b0] and return root (NaN on failure)
    /**
     * @brief Unified scalar solver interface.
     *
     * If the solver type is SimulatedAnnealing, this minimizes |f(x)| starting at x0.
     * Otherwise, it finds a root f(x) = 0 in [a0, b0].
     *
     * @param f The function to solve or minimize.
     * @param a0 Lower bound (for root finding).
     * @param b0 Upper bound (for root finding).
     * @param x0 Initial guess (for minimization).
     * @return double The result (root or minimum location). Returns NaN on failure.
     */
    double solve_scalar(const Function1D &f,
                        double a0,
                        double b0,
                        double x0) const;

    // MINIMIZATION (Simulated Annealing)
    /**
     * @brief Minimize a 1D function using the configured solver (e.g., SA).
     *
     * @param g The objective function to minimize.
     * @param x0 Initial guess.
     * @return MinResult Optimization result.
     */
    MinResult minimize(const Objective1D &g, double x0) const; // 1D

    /**
     * @brief Minimize an N-dimensional function using the configured solver (e.g., SA).
     *
     * @param g The objective function to minimize.
     * @param x0 Initial guess vector.
     * @return MinResult Optimization result.
     */
    MinResult minimize(const ObjectiveND &g,
                       const std::vector<double> &x0) const; // ND

private:
    // helpers (implemented in solver.cpp)
    MinResult sa_minimize_nd_(const ObjectiveND &g,
                            const std::vector<double> &x0) const;

private:
    CombinedSolverOptions opts_;
    Bracketer bracketer_;
};

} // namespace qss
