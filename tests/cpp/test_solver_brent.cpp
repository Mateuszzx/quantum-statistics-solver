#include <gtest/gtest.h>
#include "qss/solver.hpp"
#include <cmath>

TEST(BrentRoot, SimpleQuadratic) {
    // Solve f(x) = x^2 - 2 = 0 => x = sqrt(2)
    qss::CombinedSolverOptions opts;
    opts.solver = qss::SolverType::Brent;
    opts.root.abs_tol = 1e-12;
    opts.root.rel_tol = 1e-9;
    opts.root.max_iters = 200;
    opts.bracketer_type = qss::BracketerType::Expansion; // expansion shouldn't be needed
    qss::OptimizationSolver solver(opts);

    auto f = [](double x){ return x*x - 2.0; };
    double root = solver.solve_scalar(f, 0.0, 2.0, 0.0);
    ASSERT_TRUE(std::isfinite(root));
    ASSERT_NEAR(root, std::sqrt(2.0), 1e-9);
}

TEST(BrentRoot, SignChangeEdge) {
    // Bracket with one endpoint exactly at root
    qss::CombinedSolverOptions opts;
    opts.solver = qss::SolverType::Brent;
    opts.root.abs_tol = 1e-10;
    opts.root.rel_tol = 1e-8;
    opts.root.max_iters = 100;
    qss::OptimizationSolver solver(opts);
    auto f = [](double x){ return (x-1.25)*(x-1.25); }; // minimum at 1.25, not a root unless x==1.25 giving 0
    double root = solver.solve_scalar(f, 1.25, 2.0, 1.25);
    // Since f(a)=0 we expect to get a directly.
    ASSERT_NEAR(root, 1.25, 1e-12);
}

TEST(BrentRoot, FailureInvalidBracket) {
    // Function with same sign on interval; Brent should attempt expansion, may fail => NaN
    qss::CombinedSolverOptions opts;
    opts.solver = qss::SolverType::Brent;
    opts.root.abs_tol = 1e-12;
    opts.root.rel_tol = 1e-9;
    opts.root.max_iters = 50;
    opts.bracketer_type = qss::BracketerType::Bisection; // one-shot bisection bracket likely fails
    qss::OptimizationSolver solver(opts);
    auto f = [](double x){ return x*x + 5.0; }; // always positive
    double root = solver.solve_scalar(f, -3.0, 3.0, 0.0);
    ASSERT_TRUE(std::isnan(root));
}
