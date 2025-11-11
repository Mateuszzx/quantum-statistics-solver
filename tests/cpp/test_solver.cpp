#include <gtest/gtest.h>
#include "qss/solver.hpp"

using namespace qss;

TEST(BracketerTest, ExpansionFindsBracket)
{
    BracketerOptions opts;
    opts.expand_factor = 1.5;
    opts.max_expands = 50;
    opts.initial_step = 0.5;
    Bracketer br(BracketerType::Expansion, opts);

    // simple linear function with root at x=1.0
    Function1D f = [](double x) { return x - 1.0; };

    // start with interval that does not bracket the root
    auto B = br.bracket(f, 0.0, 0.2);
    EXPECT_TRUE(B.valid);
    EXPECT_LT(B.a, 1.0);
    EXPECT_GT(B.b, 1.0);
}

TEST(BracketerTest, BisectionMidpoint)
{
    Bracketer br(BracketerType::Bisection, BracketerOptions{});
    Function1D f = [](double x) { return x - 2.0; };
    // initial interval [1.5,2.5] already brackets root
    auto B = br.bracket(f, 1.5, 2.5);
    EXPECT_TRUE(B.valid);
    EXPECT_LE(B.a, 2.0);
    EXPECT_GE(B.b, 2.0);
}

TEST(SolverTest, BisectionConverges)
{
    CombinedSolverOptions opts;
    opts.solver = SolverType::Bisection;
    OptimizationSolver solver(opts);

    Function1D f = [](double x) { return x * x - 2.0; };
    auto R = solver.solve_root(f, 1.0, 2.0);
    EXPECT_EQ(R.status, Status::Converged);
    const double expected = std::sqrt(2.0);
    EXPECT_NEAR(R.root, expected, 1e-8);
}

TEST(SolverTest, InvalidBracket)
{
    CombinedSolverOptions opts;
    opts.solver = SolverType::Bisection;
    OptimizationSolver solver(opts);

    // no real root for this function
    Function1D f = [](double x) { return x * x + 1.0; };
    auto R = solver.solve_root(f, -1.0, 1.0);
    EXPECT_EQ(R.status, Status::InvalidBracket);
}

TEST(SolverTest, BrentConverges)
{
    CombinedSolverOptions opts;
    opts.solver = SolverType::Brent;
    OptimizationSolver solver(opts);
    // Root at x=3 for (x-3)*(x+1)=0 choose the positive root
    Function1D f = [](double x) { return (x - 3.0) * (x + 1.0); };
    auto R = solver.solve_root(f, 2.0, 4.5);
    ASSERT_EQ(R.status, Status::Converged);
    EXPECT_NEAR(R.root, 3.0, 1e-8);
}

TEST(SolverTest, SecantConverges)
{
    CombinedSolverOptions opts;
    opts.solver = SolverType::Secant;
    OptimizationSolver solver(opts);
    Function1D f = [](double x) { return std::cos(x) - x; }; // root near ~0.739085...
    auto R = solver.solve_root(f, 0.0, 1.5);
    ASSERT_TRUE(R.status == Status::Converged || R.status == Status::MaxIterations);
    EXPECT_NEAR(R.root, 0.7390851332, 1e-6);
}

TEST(SolverTest, SimulatedAnnealing1DQuadratic)
{
    CombinedSolverOptions opts;
    opts.solver = SolverType::SimulatedAnnealing;
    opts.sa.max_iters = 5000;
    opts.sa.min_iters = 200;
    opts.sa.stall_iters = 1500;
    opts.sa.T0 = 1.0;
    opts.sa.alpha = 0.995;
    OptimizationSolver solver(opts);
    Objective1D g = [](double x) { return (x - 5.0) * (x - 5.0); };
    auto MR = solver.minimize(g, 0.0);
    ASSERT_EQ(MR.status, Status::Converged);
    ASSERT_FALSE(MR.x_best.empty());
    EXPECT_NEAR(MR.x_best[0], 5.0, 0.05) << "SA should find minimum near 5";
}

TEST(SolverTest, SimulatedAnnealingBoundsReflect)
{
    CombinedSolverOptions opts;
    opts.solver = SolverType::SimulatedAnnealing;
    opts.sa.max_iters = 2000;
    opts.sa.min_iters = 100;
    opts.sa.T0 = 2.0;
    opts.sa.alpha = 0.98;
    opts.sa.bounds = {{-1.0, 1.0}}; // 1D bounds
    opts.sa.reflect_bounds = true;
    OptimizationSolver solver(opts);
    Objective1D g = [](double x) { return (x - 0.3) * (x - 0.3); };
    auto MR = solver.minimize(g, 0.95); // near upper bound
    ASSERT_EQ(MR.status, Status::Converged);
    ASSERT_FALSE(MR.x_best.empty());
    EXPECT_GE(MR.x_best[0], -1.0);
    EXPECT_LE(MR.x_best[0], 1.0);
    EXPECT_NEAR(MR.x_best[0], 0.3, 0.05);
}

