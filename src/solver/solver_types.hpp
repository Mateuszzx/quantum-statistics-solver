#pragma once
#include <functional>
#include <limits>
#include <vector>

namespace qss
{

enum class Status
{
    Converged,
    MaxIterations,
    InvalidBracket,
    DomainError,
    NotBracketed,
    Unsupported
};

// ---------- aliases ----------
using Function1D = std::function<double(double)>;  // f(x)=0 (root)
using Objective1D = std::function<double(double)>; // min g(x) (1D)
using ObjectiveND =
    std::function<double(const std::vector<double> &)>; // min g(x) (ND)

// MOVE Bracket definition here so it's complete in RootResult
struct Bracket
{
    double a{}, b{}, fa{}, fb{};
    bool valid = false; // fa*fb <= 0
};

// ---------- options ----------
struct SolverOptions
{
    double abs_tol = 1e-12;
    double rel_tol = 1e-9;
    int max_iters = 200;
};

// ---------- results ----------
struct RootResult
{
    Status status = Status::NotBracketed;
    double root = std::numeric_limits<double>::quiet_NaN();
    double froot = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    Bracket final_bracket{};
};

struct MinResult
{
    Status status = Status::Unsupported;
    std::vector<double> x_best;
    double f_best = std::numeric_limits<double>::infinity();
    int iterations = 0; // loop count (may differ from evaluations)
    int accepts = 0;
    int rejects = 0;
    double T_final = 0.0;
};

} // namespace qss
