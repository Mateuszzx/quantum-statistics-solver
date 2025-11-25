#include "qss/solver.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <iostream>


namespace qss
{

// helpers
static inline bool close_enough(double x, double y, double atol, double rtol)
{
    return std::abs(x - y) <=
           (atol + rtol * std::max(std::abs(x), std::abs(y)));
}

// -------- root solving --------
RootResult OptimizationSolver::solve_root(const Function1D &f,
                                          double a0,
                                          double b0) const
{
    RootResult R{};
    if (opts_.solver == SolverType::SimulatedAnnealing)
    {
        R.status = Status::Unsupported;
        return R;
    }

    double a = a0, b = b0;
    if (!(std::isfinite(a) && std::isfinite(b)))
    {
        R.status = Status::DomainError;
        return R;
    }
    if (a > b)
        std::swap(a, b);

    double fa = f(a), fb = f(b);
    if (!std::isfinite(fa) || !std::isfinite(fb))
    {
        R.status = Status::DomainError;
        return R;
    }

    // Ensure bracket if required. For Brent/Secant we will try expansion; for Bisection we must succeed.
    Bracket B{a, b, fa, fb, fa * fb <= 0.0};
    if (!B.valid)
    {
        B = bracketer_.bracket(f, a, b);
    }
    // If still invalid and neither endpoint is a root -> fail early for root solvers
    if (!B.valid && (opts_.solver != SolverType::SimulatedAnnealing))
    {
        R.status = Status::InvalidBracket;
        return R;
    }
    a = B.a;
    b = B.b;
    fa = B.fa;
    fb = B.fb;

    const double atol = opts_.root.abs_tol;
    const double rtol = opts_.root.rel_tol;
    const int maxI = std::max(1, opts_.root.max_iters);

    if (opts_.solver == SolverType::Bisection)
    {
        for (int i = 0; i < maxI; ++i)
        {
            R.iterations = i + 1;
            double m = 0.5 * (a + b);
            double fm = f(m);
            if (!std::isfinite(fm))
            {
                R.status = Status::DomainError;
                return R;
            }
            if (std::abs(fm) <= atol || close_enough(a, b, atol, rtol))
            {
                R.root = m;
                R.froot = fm;
                bool v = (fa * fb <= 0.0) || (std::abs(fa) <= atol) || (std::abs(fb) <= atol);
                R.final_bracket = {a, b, fa, fb, v};
                R.status = Status::Converged;
                return R;
            }
            if (fa * fm <= 0.0)
            {
                b = m;
                fb = fm;
            }
            else
            {
                a = m;
                fa = fm;
            }
        }
        R.root = 0.5 * (a + b);
        R.froot = f(R.root);
        {
            bool v = (fa * fb <= 0.0) || (std::abs(fa) <= atol) || (std::abs(fb) <= atol);
            R.final_bracket = {a, b, fa, fb, v};
        }
        R.status = Status::MaxIterations;
        return R;
    }

    // Canonical Brent implementation
    if (opts_.solver == SolverType::Brent || opts_.solver == SolverType::Secant)
    {
        double a = B.a, b = B.b;
        double fa = B.fa, fb = B.fb;
        if (!(B.valid && fa * fb <= 0.0))
        {
            R.status = Status::InvalidBracket;
            return R;
        }
        // Ensure |fa| >= |fb|
        if (std::abs(fa) < std::abs(fb)) { std::swap(a,b); std::swap(fa,fb); }
        double c = a, fc = fa;
        bool mflag = true; // last step was bisection
        double s = b, d = 0.0;
        for (int iter = 0; iter < maxI; ++iter)
        {
            R.iterations = iter + 1;
            // Convergence check
            if (std::abs(fb) <= atol || close_enough(a,b,atol,rtol))
            {
                bool v = (fa * fb <= 0.0) || (std::abs(fa) <= atol) || (std::abs(fb) <= atol);
                R.root = b; R.froot = fb; R.final_bracket = {a,b,fa,fb,v}; R.status = Status::Converged; return R;
            }
            if (fa != fc && fb != fc && opts_.solver == SolverType::Brent)
            {
                // Inverse quadratic interpolation
                s = ( a*fb*fc / ((fa - fb)*(fa - fc)) )
                  + ( b*fa*fc / ((fb - fa)*(fb - fc)) )
                  + ( c*fa*fb / ((fc - fa)*(fc - fb)) );
            }
            else
            {
                // Secant step
                s = b - fb * (b - a) / (fb - fa);
            }
            double cond1 = !((s > std::min(a,c) && s < std::max(a,c)));
            double cond2 = mflag && std::abs(s - b) >= std::abs(b - c) / 2.0;
            double cond3 = !mflag && std::abs(s - b) >= std::abs(c - d) / 2.0;
            double cond4 = mflag && std::abs(b - c) <= atol + rtol * std::abs(b);
            double cond5 = !mflag && std::abs(c - d) <= atol + rtol * std::abs(b);
            if (cond1 || cond2 || cond3 || cond4 || cond5)
            {
                s = 0.5 * (a + b);
                mflag = true;
            }
            else
            {
                mflag = false;
            }
            double fs = f(s);
            if (!std::isfinite(fs)) { R.status = Status::DomainError; return R; }
            d = c; c = b; fc = fb;
            if (fa * fs < 0)
            {
                b = s; fb = fs;
            }
            else
            {
                a = s; fa = fs;
            }
            if (std::abs(fa) < std::abs(fb)) { std::swap(a,b); std::swap(fa,fb); }
        }
        {
            bool v = (fa * fb <= 0.0) || (std::abs(fa) <= atol) || (std::abs(fb) <= atol);
            R.root = b; R.froot = fb; R.final_bracket = {a,b,fa,fb,v};

            R.status = Status::MaxIterations; return R;
        }
    }

    R.status = Status::Unsupported;
    return R;
}

// -------- unified scalar dispatch --------
double OptimizationSolver::solve_scalar(const Function1D &f,
                                        double a0,
                                        double b0,
                                        double x0) const
{
    if (opts_.solver == SolverType::SimulatedAnnealing)
    {
        auto g = [&](double x) { return std::abs(f(x)); };
        auto mr = minimize(g, x0);
        return (mr.x_best.empty() ? std::numeric_limits<double>::quiet_NaN()
                                  : mr.x_best[0]);
    }
    auto rr = solve_root(f, a0, b0);
    // Instrumentation: guard iteration sanity (env controlled)
    const char *dbg = std::getenv("QSS_DEBUG_SOLVER");
    if (dbg && (rr.iterations < 0 || rr.iterations > opts_.root.max_iters))
    {
        fprintf(stderr,
                "[qss] anomaly: iterations=%d max=%d status=%d root=%g froot=%g\n",
                rr.iterations,
                opts_.root.max_iters,
                (int)rr.status,
                rr.root,
                rr.froot);
    }
    return rr.root;
}

// -------- SA helpers used by public API --------
// Removed unused helpers sa_next_T_ and reflect_ (logic lives solely in SA impl)

// -------- public minimize routes (call the one SA impl) --------
MinResult OptimizationSolver::sa_minimize_nd_(
    const ObjectiveND &g,
    const std::vector<double> &x0) const
{
    return sa_minimize_impl(g, x0, opts_.sa);
}

MinResult OptimizationSolver::minimize(const Objective1D &g, double x0) const
{
    if (opts_.solver != SolverType::SimulatedAnnealing)
    {
        // For now only SA is supported in minimize()
        return MinResult{Status::Unsupported, std::vector<double>{},
                         std::numeric_limits<double>::infinity(), 0, 0, 0, 0.0};
    }
    return sa_minimize_nd_(
        [&](const std::vector<double> &v) { return g(v[0]); },
        std::vector<double>{x0});
}

MinResult OptimizationSolver::minimize(const ObjectiveND &g,
                                       const std::vector<double> &x0) const
{
    if (opts_.solver != SolverType::SimulatedAnnealing)
    {
        return MinResult{Status::Unsupported, std::vector<double>{},
                         std::numeric_limits<double>::infinity(), 0, 0, 0, 0.0};
    }
    return sa_minimize_nd_(g, x0);
}

} // namespace qss
