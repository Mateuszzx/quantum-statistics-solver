#include "qss/solver.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

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

    // Ensure bracket if required
    Bracket B{a, b, fa, fb, fa * fb <= 0.0};
    if (!B.valid && (opts_.solver == SolverType::Bisection ||
                     opts_.solver == SolverType::Brent))
    {
        B = bracketer_.bracket(f, a, b);
        if (!B.valid)
        {
            R.status = Status::InvalidBracket;
            return R;
        }
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
                R.final_bracket = {a, b, fa, fb, true};
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
        R.final_bracket = {a, b, fa, fb, true};
        R.status = Status::MaxIterations;
        return R;
    }

    // Brent-Dekker method
    if (opts_.solver == SolverType::Brent || opts_.solver == SolverType::Secant)
    {
        // For pure Secant, we still start from [a,b]; if not bracketed, we already expanded above.
        double c = a, fc = fa;
        double d = b - a, e = d;
        double x = b, w = b, v = b;
        double fx = fb, fw = fb, fv = fb;

        for (int i = 0; i < maxI; ++i)
        {
            R.iterations = i + 1;

            // Keep c as the previous best opposite sign point to x
            if ((fx > 0 && fc > 0) || (fx < 0 && fc < 0))
            {
                c = a;
                fc = fa;
                d = e = x - a;
            }
            if (std::abs(fc) < std::abs(fx))
            {
                a = x;
                fa = fx;
                x = c;
                fx = fc;
                c = a;
                fc = fa;
            }

            // Convergence checks
            const double xm = 0.5 * (c - x);
            const double tol1 = atol + rtol * std::abs(x);
            if (std::abs(xm) <= 2 * tol1 || std::abs(fx) <= atol)
            {
                R.root = x;
                R.froot = fx;
                R.final_bracket = {std::min(x, c),
                                   std::max(x, c),
                                   fx,
                                   fc,
                                   true};
                R.status = Status::Converged;
                return R;
            }

            // Attempt inverse quadratic interpolation or secant
            double p = 0.0, q = 1.0, r = 1.0;
            if (std::abs(e) > tol1 && std::abs(fa - fx) > 0 &&
                std::abs(fc - fx) > 0)
            {
                // Inverse quadratic interpolation using points (a,fa), (x,fx), (c,fc)
                double s = fx / fa;
                double t = fx / fc;
                double u = fx / fw; // fallback if needed
                p = s * ((c - x) * t * (t - u) - (a - x) * (1 - u));
                q = (t - 1) * (u - 1) * (s - 1);
            }
            else
            {
                // Secant step
                p = (c - x) * fx;
                q = fc - fx;
            }

            if (q > 0)
                p = -p;
            else
                q = -q;

            double min1 = 3 * xm * q - std::abs(tol1 * q);
            double min2 = std::abs(e * q);
            bool accept_interp = (std::abs(p) < std::min(min1, min2));

            if (!accept_interp)
            {
                e = d = xm;
            }
            else
            {
                e = d;
                d = p / q;
            }

            double u_new =
                (std::abs(d) > tol1) ? x + d : x + (d > 0 ? +tol1 : -tol1);
            double fu = f(u_new);
            if (!std::isfinite(fu))
            {
                R.status = Status::DomainError;
                return R;
            }

            // Update points
            a = ((fx * fu) < 0 ? x : c);
            if ((fx * fu) < 0)
            {
                c = x;
                fc = fx;
            }
            else
            {
                fc = fu;
            }
            w = x;
            fw = fx;
            x = u_new;
            fx = fu;
            fa = f(a);
        }

        R.root = x;
        R.froot = fx;
        R.final_bracket = {std::min(x, c), std::max(x, c), fx, fc, true};
        R.status = Status::MaxIterations;
        return R;
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
    return rr.root;
}

// -------- SA helpers used by public API --------
double OptimizationSolver::sa_next_T_(std::size_t k) const
{
    switch (opts_.sa.sched)
    {
    case CoolingSchedule::Exponential:
        return opts_.sa.T0 * std::pow(opts_.sa.alpha, (double)k);
    case CoolingSchedule::Linear:
    {
        if (opts_.sa.max_iters == 0)
            return 0.0;
        double frac = 1.0 - (double)k / (double)opts_.sa.max_iters;
        return std::max(0.0, opts_.sa.T0 * frac);
    }
    case CoolingSchedule::Logarithmic:
        return opts_.sa.T0 / std::log(2.0 + (double)k);
    }
    return 0.0;
}

double OptimizationSolver::reflect_(double x, double lo, double hi)
{
    if (!(std::isfinite(lo) && std::isfinite(hi)) || lo >= hi)
        return x;
    const double L = hi - lo;
    double y = std::fmod(x - lo, 2 * L);
    if (y < 0)
        y += 2 * L;
    return (y <= L) ? (lo + y) : (hi - (y - L));
}

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
        return MinResult{.status = Status::Unsupported};
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
        return MinResult{.status = Status::Unsupported};
    }
    return sa_minimize_nd_(g, x0);
}

} // namespace qss
