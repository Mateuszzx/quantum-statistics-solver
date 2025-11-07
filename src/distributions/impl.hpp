#pragma once
#include "qss/types.hpp"

namespace qss::detail
{

// ---- scalar kernels (private) ----
double fd_distribution(double e_k, double mu, double beta) noexcept;
double flb_distribution(double e_k, double mu, double beta) noexcept;
double nflb_distribution(double e_k, double mu, double beta) noexcept;
double nflb_i_distribution(double e_k, double mu, double beta) noexcept;

// ---- buffer kernels (private; fill preallocated out) ----
void fd_distribution_into(const double *ek,
                          std::size_t size,
                          double mu,
                          double beta,
                          double *out);
void flb_distribution_into(const double *ek,
                           std::size_t size,
                           double mu,
                           double beta,
                           double *out);
void nflb_distribution_into(const double *ek,
                            std::size_t size,
                            double mu,
                            double beta,
                            double *out);
void nflb_i_distribution_into(const double *ek,
                              std::size_t size,
                              double mu,
                              double beta,
                              double *out);

} // namespace qss::detail
