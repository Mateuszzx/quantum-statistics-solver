#pragma once
#include <functional>
#include <limits>
#include <vector>

namespace qss
{

/**
 * @brief Status codes for solver operations.
 */
enum class Status
{
    Converged,      ///< The solver successfully converged to a solution.
    MaxIterations,  ///< The maximum number of iterations was reached.
    InvalidBracket, ///< The provided or computed bracket is invalid (does not contain a root).
    DomainError,    ///< A domain error occurred (e.g., NaN or Inf values encountered).
    NotBracketed,   ///< The root is not bracketed within the interval.
    Unsupported     ///< The requested operation or solver type is not supported.
};

// ---------- aliases ----------
using Function1D = std::function<double(double)>;  ///< 1D function signature f(x).
using Objective1D = std::function<double(double)>; ///< 1D objective function signature min g(x).
using ObjectiveND =
    std::function<double(const std::vector<double> &)>; ///< N-dimensional objective function signature min g(x).

// MOVE Bracket definition here so it's complete in RootResult
/**
 * @brief Represents a bracketing interval [a, b] for root finding.
 */
struct Bracket
{
    double a{}, b{};   ///< The endpoints of the interval.
    double fa{}, fb{}; ///< The function values at the endpoints.
    bool valid = false; ///< True if the bracket is valid (i.e., fa * fb <= 0).
};

// ---------- options ----------
/**
 * @brief General options for root-finding solvers.
 */
struct SolverOptions
{
    double abs_tol = 1e-12; ///< Absolute tolerance for convergence.
    double rel_tol = 1e-9;  ///< Relative tolerance for convergence.
    int max_iters = 200;    ///< Maximum number of iterations allowed.
};

// ---------- results ----------
/**
 * @brief Result structure for root-finding operations.
 */
struct RootResult
{
    Status status = Status::NotBracketed; ///< The final status of the operation.
    double root = std::numeric_limits<double>::quiet_NaN(); ///< The found root, or NaN if failed.
    double froot = std::numeric_limits<double>::quiet_NaN(); ///< The function value at the root.
    int iterations = 0; ///< The number of iterations performed.
    Bracket final_bracket{}; ///< The final bracket state.
};

/**
 * @brief Result structure for minimization operations.
 */
struct MinResult
{
    Status status = Status::Unsupported; ///< The final status of the operation.
    std::vector<double> x_best; ///< The best coordinates found.
    double f_best = std::numeric_limits<double>::infinity(); ///< The minimum function value found.
    int iterations = 0; ///< The number of iterations performed.
    int accepts = 0; ///< Number of accepted moves (for SA).
    int rejects = 0; ///< Number of rejected moves (for SA).
    double T_final = 0.0; ///< Final temperature (for SA).
};

} // namespace qss
