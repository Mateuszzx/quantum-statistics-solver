#include <cmath>
#include <stdexcept>
#include <limits>
#include "qss/solver.hpp"

namespace qss {

class BrentSolver final : public IRootSolver {
public:
  double solve(const std::function<double(double)>& f,
               double a, double b,
               const RootOptions& opt) const override {
    double fa = f(a), fb = f(b);
    if (!(fa * fb <= 0.0) || !std::isfinite(fa) || !std::isfinite(fb))
      throw std::invalid_argument("Brent: invalid bracket");

    double c = a, fc = fa;
    double d = b - a, e = d;

    for (int iter = 0; iter < opt.max_iters; ++iter) {
      if (std::fabs(fc) < std::fabs(fb)) {
        std::swap(a, b); std::swap(b, c);
        std::swap(fa, fb); std::swap(fb, fc);
      }
      const double tol1 = 2.0 * std::numeric_limits<double>::epsilon() * std::fabs(b) + 0.5 * opt.tol;
      const double xm   = 0.5 * (c - b);
      if (std::fabs(xm) <= tol1 || fb == 0.0) return b;

      if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) {
        double s = fb / fa;
        double p, q;
        if (a == c) {            // secant
          p = 2.0 * xm * s;
          q = 1.0 - s;
        } else {                 // inverse quadratic
          double q1 = fa / fc;
          double r  = fb / fc;
          p = s * (2.0 * xm * q1 * (q1 - r) - (b - a) * (r - 1.0));
          q = (q1 - 1.0) * (r - 1.0) * (s - 1.0);
        }
        if (p > 0) q = -q; else p = -p;

        const double min1 = 3.0 * xm * q - std::fabs(tol1 * q);
        const double min2 = std::fabs(e * q);
        if (2.0 * p < (min1 < min2 ? min1 : min2)) { e = d; d = p / q; }
        else { d = xm; e = d; }
      } else {
        d = xm; e = d;
      }

      a = b; fa = fb;
      b += (std::fabs(d) > tol1) ? d : (xm >= 0 ? tol1 : -tol1);
      fb = f(b);
      if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; e = d = b - a; }
    }
    return b; // best effort (or throw if you prefer)
  }
};

std::unique_ptr<IRootSolver> make_solver(SolverType t) {
  switch (t) {
    case SolverType::Brent: return std::make_unique<BrentSolver>();
    default:                return std::make_unique<BrentSolver>();
  }
}

} // namespace qss
