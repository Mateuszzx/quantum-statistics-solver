#pragma once
#include "qss/solver.hpp" // for enums

namespace qss
{

struct SolverParams
{
    SolverType solver = SolverType::Brent; // choose implementation
    BracketerType bracketer = BracketerType::Expansion;
    double tol = 1e-6; // RootOptions
    int max_iters = 300;
    double expand_factor = 1.5; // BracketerOptions
    int max_expands = 200;
    double a0 = -50.0; // initial window
    double b0 = 10.0;
};

} // namespace qss
