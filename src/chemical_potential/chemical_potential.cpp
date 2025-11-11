// Simplified chemical potential solving implementation
#include "qss/chemical_potential.hpp"
#include "chem_helpers.hpp" // distribution_total(E, mu, beta, dist, delta)
#include "qss/math_utils.hpp"
#include "qss/solver.hpp"
#include <algorithm>
#include <limits>
#include <optional>
#include <stdexcept>

namespace qss
{

namespace
{
    inline double zeeman_delta(double H_a_T, double g, double muB_eV_per_T)
    {
        return 0.5 * g * muB_eV_per_T * H_a_T; // same units as E
    }

    // Scan [lo,hi] in "steps" increments to find a bracket where f changes sign
    static std::optional<std::pair<double, double>> scan_bracket(
        const std::function<double(double)> &f, double lo, double hi, int steps)
    {
        if (!(std::isfinite(lo) && std::isfinite(hi)) || !(lo < hi) || steps < 2)
            return std::nullopt;
        double prev_x = lo;
        double prev_f = f(prev_x);
        if (!std::isfinite(prev_f))
            return std::nullopt;
        for (int i = 1; i <= steps; ++i)
        {
            double x = lo + (hi - lo) * (double)i / (double)steps;
            double fx = f(x);
            if (!std::isfinite(fx))
                continue;
            if (prev_f == 0.0)
                return std::make_pair(prev_x, prev_x);
            if (fx == 0.0)
                return std::make_pair(x, x);
            if (prev_f * fx < 0.0)
                return std::make_pair(prev_x, x);
            prev_x = x;
            prev_f = fx;
        }
        return std::nullopt;
    }

    // Unified internal solver returning diagnostics (callers can discard fields)
    static ChemicalPotentialResult solve_mu(
        DistType dist,
        double beta,
        double n_assumed,
        const double N,
        const Mat &E,
        double delta,
        const CombinedSolverOptions &opts,
        double a0 = -5.0,
        double b0 = 5.0,
        double x0 = 0.0)
    {
        if (E.empty() || E[0].empty())
            throw std::invalid_argument("chemical_potential_B: empty E");
        if (N != E.size() || N != E[0].size())
            throw std::invalid_argument("chemical_potential_B: kx/ky vs E shape mismatch");

        auto f = [&](double mu) {
            return qss::math::mean_abs_2d(qss::chem::distribution_total(E, mu, beta, dist, delta), N) - n_assumed;
        };

        qss::OptimizationSolver solver(opts);
        ChemicalPotentialResult R{};


        if (opts.solver == qss::SolverType::SimulatedAnnealing)
        {
            auto g = [&](double x) { return std::abs(f(x)); };
            CombinedSolverOptions sa_opts = opts;
            if (sa_opts.sa.bounds.empty())
            {
                sa_opts.sa.bounds = {std::make_pair(a0, b0)};
            }

            qss::OptimizationSolver sa_solver(sa_opts);
            MinResult mr = sa_solver.minimize(g, x0);

            if (!mr.x_best.empty())
            {
                R.mu = mr.x_best[0];
            }

            R.iterations = static_cast<int>(mr.iterations);
            R.converged = (mr.status == qss::Status::Converged);
        }
        else
        {
            RootResult rr = solver.solve_root(f, a0, b0);
            if (rr.status != qss::Status::Converged)
            {
                // coarse scan for a valid sign change
                if (auto B = scan_bracket(f, a0, b0, 96))
                {
                    rr = solver.solve_root(f, B->first, B->second);
                }
            }
            if (rr.status == qss::Status::Converged)
            {
                R.mu = rr.root;
                R.iterations = static_cast<int>(rr.iterations);
                R.converged = true;
            }
        }
        if (std::isfinite(R.mu))
            R.error = std::abs(f(R.mu));
        return R;
    }
} // namespace

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
    double a0,
    double b0,
    double x0)
{
    const double delta = 0.5 * g * muB_eV_per_T * H_a_T;
    ChemicalPotentialResult R = solve_mu(
        dist, beta, n_assumed, N, E, delta, solver_opts, a0, b0, x0);

    if (std::isfinite(R.mu)) {
        const auto F = qss::chem::distribution_total(E, R.mu, beta, dist, delta);
        R.error = std::abs(qss::math::mean_abs_2d(F, N) - n_assumed);
    }
    return R;
}

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
    double a0,
    double b0,
    double x0)
{
    std::vector<ChemicalPotentialResult> out;
    out.reserve(betas.size());
    for (double b : betas) {
        out.push_back(chemical_potential_B(
            dist, b, n_assumed, N, E, H_a_T, g, muB_eV_per_T,
            solver_opts, a0, b0, x0));
    }
    return out;
}

// Multi-n
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
    double a0,
    double b0,
    double x0)
{
    std::vector<ChemicalPotentialResult> out;
    out.reserve(n_assumed_values.size());
    for (double n_assumed : n_assumed_values) {
        out.push_back(chemical_potential_B(
            dist,
            beta,
            n_assumed,
            N,
            E,
            H_a_T,
            g,
            muB_eV_per_T,
            solver_opts,
            a0,
            b0,
            x0));
    }
    return out;
}

} // namespace qss
