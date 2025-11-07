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

    if (type_ == BracketerType::Expansion)
    {
        double ef = std::max(1.0001, options().expand_factor);
        int k = std::max(0, options().max_expands);
        double step =
            (options().initial_step > 0 ? options().initial_step : 1.0);

        double left = a, right = b;
        double fl = fa, fr = fb;

        while (k-- > 0)
        {
            // expand towards the side with smaller |f|
            if (std::abs(fl) < std::abs(fr))
            {
                left -= step;
                fl = f(left);
            }
            else
            {
                right += step;
                fr = f(right);
            }
            if (!std::isfinite(fl) || !std::isfinite(fr))
                break;
            if (fl * fr <= 0.0)
            {
                B = {left, right, fl, fr, true};
                return B;
            }
            step *= ef;
        }
    }
    return B; // invalid
}

} // namespace qss
