/**
 * @file distributions.hpp
 * @brief Defines statistical distribution functions for quantum systems.
 *
 * This file contains the definitions for various statistical distributions used
 * in quantum statistics calculations, including Fermi-Dirac and non-Fermi-liquid
 * variants. It provides a unified interface via the `distribution` function
 * overloads.
 */

#pragma once
#include "types.hpp"
#include <optional>
#include <tuple>

namespace qss
{

/**
 * @brief Enumeration of supported statistical distribution types.
 */
enum class DistType
{
    FD,     ///< Fermi-Dirac distribution
    FLB,    ///< Fermi-Liquid-Like distribution
    NFL,    ///< Non-Fermi-Liquid distribution
    NFLB,   ///< Non-Fermi-Liquid with Boltzmann distribution
    SSL,   ///< Spin-Split Liquid distribution
    NFLB_i  ///< Implicit Non-Fermi-Liquid with Boltzmann distribution
};


/**
 * @brief Calculate the distribution function value for a single energy state.
 *
 * @param e_k Energy of the state.
 * @param mu Chemical potential.
 * @param beta Inverse temperature (1/kT).
 * @param delta Zeeman shift.
 * @param dist Distribution type.
 * @return std::optional<double> Distribution value, or std::nullopt if dist is invalid.
 */
std::optional<double> distribution(double e_k,
                                   double mu,
                                   double beta,
                                   double delta,
                                   int spin_shift,
                                   DistType dist) noexcept;

/**
 * @brief Calculate the distribution function for a vector of energy states.
 *
 * @param e_k Vector of energy states.
 * @param mu Chemical potential.
 * @param beta Inverse temperature (1/kT).
 * @param delta Zeeman shift.
 * @param dist Distribution type.
 * @return std::optional<Vec> Vector of distribution values, or std::nullopt if dist is invalid.
 */
std::optional<Vec> distribution(const Vec &e_k,
                                double mu,
                                double beta,
                                double delta,
                                int spin_shift,
                                DistType dist);

/**
 * @brief Calculate the distribution function for a 2D matrix of energy states.
 *
 * @param e_k Matrix of energy states.
 * @param mu Chemical potential.
 * @param beta Inverse temperature (1/kT).
 * @param delta Zeeman shift.
 * @param dist Distribution type.
 * @return std::optional<Mat> Matrix of distribution values, or std::nullopt if dist is invalid.
 */
std::optional<Mat> distribution(const Mat &e_k,
                                double mu,
                                double beta,
                                double delta,
                                int spin_shift,
                                DistType dist);

} // namespace qss

namespace qss::detail
{

// ---- scalar kernels (private) ----

/**
 * @brief Compute Fermi-Dirac distribution for a single state.
 * @internal
 */
double fd_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept;

/**
 * @brief Compute Fermi-Liquid-Like distribution for a single state.
 * @internal
 */
double flb_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept;

/**
 * @brief Compute Non-Fermi-Liquid distribution for a single state.
 * @internal
 */
double nfl_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept;

/**
 * @brief Compute Non-Fermi-Liquid with Boltzmann distribution for a single state.
 * @internal
 */
double nflb_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept;

/**
 * @brief Compute Spin-Split Liquid for a single state.
 * @internal
 */
double ssl_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept;

/**
 * @brief Compute Implicit Non-Fermi-Liquid with Boltzmann distribution for a single state.
 * @internal
 */
double nflb_i_distribution(double e_k, double mu, double beta, double delta, int spin_shift) noexcept;

} // namespace qss::detail
