#pragma once
#include <functional>
#include <memory>
#include <utility>

namespace qss {

// -------- Bracketer interface --------
struct Bracket {
  double a{};
  double b{};
  bool   ok{false};
};

struct BracketerOptions {
  double expand_factor = 1.5;
  int    max_expands   = 200;
};

class IBracketer {
public:
  virtual ~IBracketer() = default;
  // Try to find [a,b] with f(a)*f(b) < 0. Returns {a,b,true} on success.
  virtual Bracket bracket(const std::function<double(double)>& f,
                          double a0, double b0,
                          const BracketerOptions& opt) const = 0;
};

// -------- Root solver interface --------
struct RootOptions {
  double tol       = 1e-6;
  int    max_iters = 300;
};

class IRootSolver {
public:
  virtual ~IRootSolver() = default;
  // Solve f(x)=0 on [a,b] with f(a)*f(b) < 0. Returns root.
  virtual double solve(const std::function<double(double)>& f,
                       double a, double b,
                       const RootOptions& opt) const = 0;
};

// -------- Factories (optional) --------
enum class SolverType { Brent /*, Secant, Bisection, ...*/ };
enum class BracketerType { Expansion /*, ExpansionContraction, ...*/ };

std::unique_ptr<IRootSolver>   make_solver(SolverType);
std::unique_ptr<IBracketer>    make_bracketer(BracketerType);

} // namespace qss
