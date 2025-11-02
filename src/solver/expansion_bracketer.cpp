#include <cmath>
#include <algorithm>
#include "qss/solver.hpp"

namespace qss {

class ExpansionBracketer final : public IBracketer {
public: 
  Bracket bracket(const std::function<double(double)>& f,
                  double a0, double b0,
                  const BracketerOptions& opt) const override {
    if (a0 >= b0) std::swap(a0, b0);
    double fa = f(a0), fb = f(b0);
    if (!std::isfinite(fa) || !std::isfinite(fb)) return {};
    if (fa * fb < 0.0) return {a0, b0, true};

    const double mid0 = 0.5 * (a0 + b0);
    double radius     = 0.5 * (b0 - a0);

    for (int k = 0; k < opt.max_expands; ++k) {
      radius *= opt.expand_factor;
      const double a = mid0 - radius;
      const double b = mid0 + radius;
      const double fa2 = f(a), fb2 = f(b);
      if (std::isfinite(fa2) && std::isfinite(fb2) && fa2 * fb2 < 0.0)
        return {a, b, true};
    }
    return {}; // ok=false
  }
};

std::unique_ptr<IBracketer> make_bracketer(BracketerType t) {
  switch (t) {
    case BracketerType::Expansion: return std::make_unique<ExpansionBracketer>();
    default:                       return std::make_unique<ExpansionBracketer>();
  }
}

} // namespace qss
