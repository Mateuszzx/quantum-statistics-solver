#pragma once
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "../src/solver/bracket.hpp"
#include "../src/solver/solver_sa.hpp"

namespace qss
{

// ---------- enums ----------
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

// Unified options: contains both bracketer options and per-solver child options
struct CombinedSolverOptions
{
    SolverType solver = SolverType::Brent;
    BracketerType bracketer_type = BracketerType::Expansion;
    BracketerOptions bracketer_opts{};
    RootSolverOptions root{};       // child options for root solvers
    SimulatedAnnealingOptions sa{}; // child options for SA
};

// ---------- OptimizationSolver (public interface) ----------
class OptimizationSolver
{
public:
    explicit OptimizationSolver(const CombinedSolverOptions &opts = {})
        : opts_(opts), bracketer_(opts_.bracketer_type, opts_.bracketer_opts)
    {
    }

    // configuration
    void set_options(const CombinedSolverOptions &o)
    {
        opts_ = o;
        bracketer_ = Bracketer(opts_.bracketer_type, opts_.bracketer_opts);
    }

    void set_solver_type(SolverType t) noexcept
    {
        opts_.solver = t;
    }
    void set_bracketer(BracketerType t, const BracketerOptions &o = {}) noexcept
    {
        opts_.bracketer_type = t;
        opts_.bracketer_opts = o;
        bracketer_ = Bracketer(t, o);
    }
    void set_root_options(const RootSolverOptions &o) noexcept
    {
        opts_.root = o;
    }
    void set_sa_options(const SimulatedAnnealingOptions &o) noexcept
    {
        opts_.sa = o;
    }

    const CombinedSolverOptions &options() const noexcept
    {
        return opts_;
    }
    const Bracketer &bracketer() const noexcept
    {
        return bracketer_;
    }

    // ROOT-FINDING (bisection/secant/brent)
    RootResult solve_root(const Function1D &f, double a0, double b0) const;

    // Unified scalar solve:
    // - if solver == SimulatedAnnealing -> minimize |f(x)| starting at x0 and return x*
    // - else -> root-find f(x)=0 on [a0,b0] and return root (NaN on failure)
    double solve_scalar(const Function1D &f,
                        double a0,
                        double b0,
                        double x0) const;

    // MINIMIZATION (Simulated Annealing)
    MinResult minimize(const Objective1D &g, double x0) const; // 1D
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
