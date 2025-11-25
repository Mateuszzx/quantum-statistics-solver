#include <gtest/gtest.h>
#include "qss/solver.hpp"
#include "qss/solver_sa.hpp"
#include <cmath>
#include <limits>

using namespace qss;

// --- Helper Functions ---
static double simple_parabola(double x) { return x * x; }
static double shifted_parabola(double x) { return (x - 2.0) * (x - 2.0); }

TEST(SimulatedAnnealingTest, 1DQuadratic)
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

TEST(SimulatedAnnealingTest, BoundsReflect)
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

TEST(SimulatedAnnealingTest, SolveScalarDispatch) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::SimulatedAnnealing);
    
    SimulatedAnnealingOptions sa_opts;
    sa_opts.T0 = 1.0;
    sa_opts.max_iters = 1000;
    sa_opts.step_sigma = 0.5;
    solver.set_sa_options(sa_opts);
    
    // Minimizing (x-2)^2, minimum at x=2
    // solve_scalar with SA performs minimization
    double result = solver.solve_scalar(shifted_parabola, 0.0, 5.0, 2.5);
    
    // Result should be close to 2.0
    EXPECT_NEAR(result, 2.0, 0.1); 
}

TEST(SimulatedAnnealingTest, MinimizeND_Simple) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::SimulatedAnnealing);
    
    SimulatedAnnealingOptions sa_opts;
    sa_opts.max_iters = 5000;
    sa_opts.T0 = 10.0;
    sa_opts.step_sigma = 2.0; // Larger steps
    sa_opts.step_decay = 0.999; // Slow decay
    solver.set_sa_options(sa_opts);
    
    // Minimize x^2 + y^2
    ObjectiveND func = [](const std::vector<double>& x) {
        return x[0]*x[0] + x[1]*x[1];
    };
    
    std::vector<double> x0 = {5.0, -5.0};
    auto result = solver.minimize(func, x0);
    
    EXPECT_EQ(result.status, Status::Converged);
    EXPECT_NEAR(result.f_best, 0.0, 0.5); // Looser tolerance
    EXPECT_NEAR(result.x_best[0], 0.0, 0.5);
    EXPECT_NEAR(result.x_best[1], 0.0, 0.5);
}

TEST(SimulatedAnnealingTest, MinimizeND_EmptyX0) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::SimulatedAnnealing);
    
    ObjectiveND func = [](const std::vector<double>& x) { return 0.0; };
    std::vector<double> x0; // Empty
    
    auto result = solver.minimize(func, x0);
    EXPECT_EQ(result.status, Status::DomainError);
}

TEST(SimulatedAnnealingTest, CoolingSchedules) {
    OptimizationSolver solver;
    solver.set_solver_type(SolverType::SimulatedAnnealing);
    
    SimulatedAnnealingOptions sa_opts;
    sa_opts.max_iters = 100;
    
    // Linear
    sa_opts.sched = CoolingSchedule::Linear;
    solver.set_sa_options(sa_opts);
    double res_lin = solver.solve_scalar(simple_parabola, -5.0, 5.0, 2.0);
    EXPECT_TRUE(std::isfinite(res_lin));

    // Logarithmic
    sa_opts.sched = CoolingSchedule::Logarithmic;
    solver.set_sa_options(sa_opts);
    double res_log = solver.solve_scalar(simple_parabola, -5.0, 5.0, 2.0);
    EXPECT_TRUE(std::isfinite(res_log));
}

TEST(SimulatedAnnealingTest, BoundsClamp) {
    // Test the non-reflecting bounds logic
    SimulatedAnnealingOptions opts;
    opts.max_iters = 200;
    opts.reflect_bounds = false; // Clamp
    opts.bounds = {{-1.0, 1.0}}; // 1D bounds
    opts.T0 = 50.0;
    opts.step_sigma = 2.0; // Large steps to hit bounds easily
    
    // Function that decreases as x goes to +infinity, so it pushes against upper bound
    ObjectiveND func = [](const std::vector<double>& x) {
        return -x[0]; 
    };
    
    std::vector<double> x0 = {0.0};
    auto result = sa_minimize_impl(func, x0, opts);
    
    // Should be clamped at 1.0
    EXPECT_NEAR(result.x_best[0], 1.0, 1e-2);
}

TEST(SimulatedAnnealingTest, InfiniteFunction) {
    SimulatedAnnealingOptions opts;
    opts.max_iters = 50;
    opts.T0 = 1.0;
    
    // Function returns infinity sometimes
    ObjectiveND func = [](const std::vector<double>& x) {
        if (std::abs(x[0]) < 0.5) return std::numeric_limits<double>::infinity();
        return x[0]*x[0];
    };
    
    std::vector<double> x0 = {2.0};
    auto result = sa_minimize_impl(func, x0, opts);
    
    // Should not crash, should find a valid minimum (likely near original if blocked)
    EXPECT_TRUE(std::isfinite(result.f_best));
}

TEST(SimulatedAnnealingTest, InitialInfinite) {
    SimulatedAnnealingOptions opts;
    opts.max_iters = 10;
    
    // Initial point is infinite
    ObjectiveND func = [](const std::vector<double>& x) {
        return std::numeric_limits<double>::infinity();
    };
    
    std::vector<double> x0 = {0.0};
    auto result = sa_minimize_impl(func, x0, opts);
    
    EXPECT_EQ(result.status, Status::DomainError);
}
