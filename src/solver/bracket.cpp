#include "bracket.hpp"
#include <algorithm>
#include <cmath>

namespace qss
{

Bracket Bracketer::bracket(const Function1D &f, double a0, double b0) const
{
    Bracket B{};
    double a = a0, b = b0;
    if (a > b)
        std::swap(a, b);
    double fa = f(a), fb = f(b);

    if (!std::isfinite(fa) || !std::isfinite(fb))
        return B;
    if (fa == 0.0)
    {
        B = {a, a, fa, fa, true};
        return B;
    }
    if (fb == 0.0)
    {
        B = {b, b, fb, fb, true};
        return B;
    }
    if (fa * fb <= 0.0)
    {
        B = {a, b, fa, fb, true};
        return B;
    }

    // Cache options to avoid repeated accessor calls and improve clarity
    const BracketerOptions opt = options();

    if (type_ == BracketerType::Expansion)
    {
        // Quick midpoint check first
        double mid = 0.5 * (a + b);
        double fm = f(mid);
        if (std::isfinite(fm))
        {
            if (fa * fm <= 0.0)
                return Bracket{a, mid, fa, fm, true};
            if (fm * fb <= 0.0)
                return Bracket{mid, b, fm, fb, true};
        }

        // Pre-scan uniformly for an early bracket (helps with flat plateaus)
        const int scan_pts = 64;
        double prev_x = a;
        double prev_f = fa;
        for (int i = 1; i <= scan_pts; ++i)
        {
            double x = a + (b - a) * (double)i / (double)scan_pts;
            double fx = f(x);
            if (!std::isfinite(fx))
                continue;
            if (prev_f * fx <= 0.0)
            {
                return Bracket{prev_x, x, prev_f, fx, true};
            }
            prev_x = x;
            prev_f = fx;
        }

        // Symmetric expansion about the original endpoints.
        double ef = std::max(1.0001, opt.expand_factor);
        int k = std::max(0, opt.max_expands);
        double step = (opt.initial_step > 0 ? opt.initial_step : 1.0);

        double left = a;
        double right = b;
        double fl = fa;
        double fr = fb;

        while (k-- > 0)
        {
            left -= step;
            right += step;
            fl = f(left);
            fr = f(right);
            if (!std::isfinite(fl) || !std::isfinite(fr))
                break;

            if (fl * fr <= 0.0)
            {
                return Bracket{left, right, fl, fr, true};
            }
            double mid2 = 0.5 * (left + right);
            double fm2 = f(mid2);
            if (std::isfinite(fm2))
            {
                if (fl * fm2 <= 0.0)
                    return Bracket{left, mid2, fl, fm2, true};
                if (fm2 * fr <= 0.0)
                    return Bracket{mid2, right, fm2, fr, true};
            }
            step *= ef;
        }
    }
    else if (type_ == BracketerType::Bisection)
    {
        // Simple bisection-based bracketing: check the midpoint for a sign change
        // This does a single subdivision and returns a valid bracket if one
        // of the halves contains a sign change. It does not iterate further.
        double mid = 0.5 * (a + b);
        double fm = f(mid);
        if (std::isfinite(fm))
        {
            if (fa * fm <= 0.0)
            {
                B = {a, mid, fa, fm, true};
                return B;
            }
            if (fm * fb <= 0.0)
            {
                B = {mid, b, fm, fb, true};
                return B;
            }
        }
        // if midpoint didn't bracket, return invalid result
    }
    return B; // invalid
}

} // namespace qss
