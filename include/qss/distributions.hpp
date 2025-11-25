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
    NFLB,   ///< Non-Fermi-Liquid distribution
    NFLB_i  ///< Implicit Non-Fermi-Liquid distribution
};


/**
 * @brief Calculate the distribution function value for a single energy state.
 *
 * @param e_k Energy of the state.
 * @param mu Chemical potential.
 * @param beta Inverse temperature (1/kT).
 * @param dist Distribution type.
 * @return std::optional<double> Distribution value, or std::nullopt if dist is invalid.
 */
std::optional<double> distribution(double e_k,
                                   double mu,
                                   double beta,
                                   DistType dist) noexcept;

/**
 * @brief Calculate the distribution function for a vector of energy states.
 *
 * @param e_k Vector of energy states.
 * @param mu Chemical potential.
 * @param beta Inverse temperature (1/kT).
 * @param dist Distribution type.
 * @return std::optional<Vec> Vector of distribution values, or std::nullopt if dist is invalid.
 */
std::optional<Vec> distribution(const Vec &e_k,
                                double mu,
                                double beta,
                                DistType dist);

/**
 * @brief Calculate the distribution function for a 2D matrix of energy states.
 *
 * @param e_k Matrix of energy states.
 * @param mu Chemical potential.
 * @param beta Inverse temperature (1/kT).
 * @param dist Distribution type.
 * @return std::optional<Mat> Matrix of distribution values, or std::nullopt if dist is invalid.
 */
std::optional<Mat> distribution(const Mat &e_k,
                                double mu,
                                double beta,
                                DistType dist);

} // namespace qss

namespace qss::detail
{

// ---- scalar kernels (private) ----

/**
 * @brief Compute Fermi-Dirac distribution for a single state.
 * @internal
 */
double fd_distribution(double e_k, double mu, double beta) noexcept;

/**
 * @brief Compute Fermi-Liquid-Like distribution for a single state.
 * @internal
 */
double flb_distribution(double e_k, double mu, double beta) noexcept;

/**
 * @brief Compute Non-Fermi-Liquid distribution for a single state.
 * @internal
 */
double nflb_distribution(double e_k, double mu, double beta) noexcept;

/**
 * @brief Compute Implicit Non-Fermi-Liquid distribution for a single state.
 * @internal
 */
double nflb_i_distribution(double e_k, double mu, double beta) noexcept;

// ---- buffer kernels (private; fill preallocated out) ----

/**
 * @brief Vectorized Fermi-Dirac distribution.
 * @param ek Pointer to energy array.
 * @param size Number of elements.
 * @param mu Chemical potential.
 * @param beta Inverse temperature.
 * @param out Pointer to output array (must be preallocated).
 * @internal
 */
void fd_distribution_into(const double *ek,
                          std::size_t size,
                          double mu,
                          double beta,
                          double *out);

/**
 * @brief Vectorized Fermi-Liquid-Like distribution.
 * @see fd_distribution_into
 * @internal
 */
void flb_distribution_into(const double *ek,
                           std::size_t size,
                           double mu,
                           double beta,
                           double *out);

/**
 * @brief Vectorized Non-Fermi-Liquid distribution.
 * @see fd_distribution_into
 * @internal
 */
void nflb_distribution_into(const double *ek,
                            std::size_t size,
                            double mu,
                            double beta,
                            double *out);

/**
 * @brief Vectorized Implicit Non-Fermi-Liquid distribution.
 * @see fd_distribution_into
 * @internal
 */
void nflb_i_distribution_into(const double *ek,
                              std::size_t size,
                              double mu,
                              double beta,
                              double *out);

} // namespace qss::detail
