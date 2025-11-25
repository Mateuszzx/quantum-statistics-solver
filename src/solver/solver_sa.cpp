#include "qss/solver_sa.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace qss
{

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

static inline double next_T(const SimulatedAnnealingOptions &o, std::size_t k)
{
    switch (o.sched)
    {
    case CoolingSchedule::Exponential:
        return o.T0 * std::pow(o.alpha, (double)k);
    case CoolingSchedule::Linear:
    {
        if (o.max_iters == 0)
            return 0.0;
        double frac = 1.0 - (double)k / (double)o.max_iters;
        return std::max(0.0, o.T0 * frac);
    }
    case CoolingSchedule::Logarithmic:
        return o.T0 / std::log(2.0 + (double)k);
    }
    return 0.0;
}

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

    R.x_best = x0;

    // bring x0 into bounds if needed
    if (!o.bounds.empty())
    {
        if (o.bounds.size() != x0.size())
            throw std::runtime_error(
                "SA: bounds size must match x0 size or be empty");
        for (std::size_t d = 0; d < x0.size(); ++d)
        {
            auto [lo, hi] = o.bounds[d];
            if (std::isfinite(lo) && R.x_best[d] < lo)
                R.x_best[d] =
                    o.reflect_bounds ? reflect_box(R.x_best[d], lo, hi) : lo;
            if (std::isfinite(hi) && R.x_best[d] > hi)
                R.x_best[d] =
                    o.reflect_bounds ? reflect_box(R.x_best[d], lo, hi) : hi;
        }
    }

    double f_best = g(R.x_best);
    if (!std::isfinite(f_best))
    {
        R.status = Status::DomainError;
        return R;
    }
    R.f_best = f_best;

    std::vector<double> x_cur = R.x_best;
    double f_cur = f_best;

    std::vector<double> x_cand(x0.size());
    std::size_t no_improve = 0;
    double T = o.T0;
    const double base_sigma = std::max(1e-16, o.step_sigma);
    double sigma = base_sigma;

    for (std::size_t k = 0; k < o.max_iters; ++k)
    {
        R.iterations = k + 1;

        // temperature & step update
        T = next_T(o, k);

        // cumulative decay across iterations: base * scale(T) * (step_decay^k)
        sigma = base_sigma;
        if (o.scale_sigma_by_T && o.T0 > 0.0)
        {
            sigma *= std::sqrt(std::max(0.0, T / o.T0));
        }
        const double decay_k = std::pow(std::max(0.0, o.step_decay), static_cast<double>(k));
        sigma *= decay_k;
        sigma = std::max(1e-16, sigma);

        if (k + 1 >= o.min_iters)
        {
            if (T < o.T_min)
                break;
            if (no_improve >= o.stall_iters)
                break;
        }

        // propose candidate
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

        double dE = f_cand - f_cur;
        bool accept =
            (dE <= 0.0) || (uni(rng) < std::exp(-dE / std::max(T, 1e-16)));

        if (accept)
        {
            x_cur = x_cand;
            f_cur = f_cand;
            ++R.accepts;
            if (f_cur < f_best)
            {
                f_best = f_cur;
                R.x_best = x_cur;
                R.f_best = f_best;
                no_improve = 0;
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
