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


// --- Helper Functions ---
static double simple_parabola(double x) { return x * x; }
static double simple_root_func(double x) { return x * x - 4.0; } // Roots at +/- 2

// --- OptimizationSolver::solve_root Coverage ---

TEST(SolverTest, SolveRoot_UnsupportedSolver) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::SimulatedAnnealing); // Not a root solver
    
    auto result = solver.solve_root(simple_root_func, 0.0, 5.0);
    EXPECT_EQ(result.status, Status::Unsupported);
}

TEST(SolverTest, SolveRoot_DomainError_InfiniteBounds) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::Bisection);
    
    double inf = std::numeric_limits<double>::infinity();
    
    auto result1 = solver.solve_root(simple_root_func, -inf, 5.0);
    EXPECT_EQ(result1.status, Status::DomainError);

    auto result2 = solver.solve_root(simple_root_func, 0.0, inf);
    EXPECT_EQ(result2.status, Status::DomainError);
}

TEST(SolverTest, SolveRoot_MaxIterations_Bisection) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::Bisection);
    
    SolverOptions root_opts;
    root_opts.max_iters = 1; // Too few to converge
    root_opts.rel_tol = 1e-15;
    root_opts.abs_tol = 1e-15;
    solver.set_root_options(root_opts);

    // Root is at 2.0. Bracket [0, 5]. Midpoint 2.5. 
    // 1 iteration won't reach tolerance.
    auto result = solver.solve_root(simple_root_func, 0.0, 5.0);
    EXPECT_EQ(result.status, Status::MaxIterations);
}

TEST(SolverTest, SolveRoot_MaxIterations_Brent) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::Brent);
    
    SolverOptions root_opts;
    root_opts.max_iters = 1; 
    root_opts.rel_tol = 1e-15;
    root_opts.abs_tol = 1e-15;
    solver.set_root_options(root_opts);

    auto result = solver.solve_root(simple_root_func, 0.0, 5.0);
    EXPECT_EQ(result.status, Status::MaxIterations);
}

TEST(SolverTest, Minimize_UnsupportedSolver) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::Bisection); // Not a minimizer
    
    auto result = solver.minimize(simple_parabola, 2.5);
    EXPECT_EQ(result.status, Status::Unsupported);
}

TEST(BracketerTest, MaxIterations) {
    BracketerOptions opts;
    opts.max_expands = 1; // Very few expands
    Bracketer bracketer(BracketerType::Expansion, opts);
    
    // Start with a small interval that doesn't bracket the root (root at 2.0)
    // [0.1, 0.2] -> f(0.1) < 0, f(0.2) < 0. Needs expansion.
    auto res = bracketer.bracket(simple_root_func, 0.1, 0.2);
    
    // It might fail to find bracket in 1 step
    // If it fails, valid should be false
    if (!res.valid) {
        SUCCEED();
    } else {
        // If it succeeded, check if it actually brackets
        EXPECT_LE(res.fa * res.fb, 0.0);
    }
}

