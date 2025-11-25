#pragma once
#include <stdexcept>
#include "qss/distributions.hpp"
#include "qss/types.hpp"

namespace qss::chem {

/**
 * @brief Computes the total occupation (sum over spins) for the given parameters.
 * 
 * @param E Energy grid 2D (ny rows x nx columns).
 * @param mu Chemical potential.
 * @param beta Inverse temperature.
 * @param dist Distribution type (FD/FLB/NFLB/...).
 * @param delta Zeeman shift (g*muB*H/2) in energy units.
 * @return Mat Total distribution matrix.
 * 
 * Principle:
 * - For delta == 0: No spin splitting -> F = f(E)
 * - For delta != 0: Spin splitting -> F = f(E - delta) + f(E + delta)
 *   (summing up and down channels point-by-point)
 */
inline Mat distribution_total(const Mat& E, double mu, double beta,
                            DistType dist, double delta) {
    // Dimensions of the energy matrix
    const std::size_t ny = E.size();
    const std::size_t nx = ny ? E[0].size() : 0;

    // Simple row consistency check (optional; O(ny))
    for (std::size_t j = 1; j < ny; ++j) {
        if (E[j].size() != nx)
            throw std::runtime_error("distribution_total: E has ragged rows");
    }

    // Output has the same shape as E
    Mat out(ny, Vec(nx));

    for (std::size_t j = 0; j < ny; ++j) {
        // First calculate distribution for row E[j] without shift
        // (useful when delta == 0)
        auto f0 = qss::distribution(E[j], mu, beta, dist);
        if (!f0)
            throw std::runtime_error("distribution() failed");

        if (delta == 0.0) {
            // No splitting: copy f(E) results to output
            out[j] = std::move(*f0);
            continue; // next row
        }

        // Spin splitting: prepare two "energy rows"
        // shifted by ±delta and calculate distributions for both channels.
        Vec row_up(nx), row_dn(nx);
        for (std::size_t i = 0; i < nx; ++i) {
            row_up[i] = E[j][i] - delta; // spin ↑
            row_dn[i] = E[j][i] + delta; // spin ↓
        }

        auto f_up = qss::distribution(row_up, mu, beta, dist);
        auto f_dn = qss::distribution(row_dn, mu, beta, dist);
        if (!f_up || !f_dn)
            throw std::runtime_error("distribution() failed (split)");

        // Sum of spin channels point-by-point
        for (std::size_t i = 0; i < nx; ++i)
            out[j][i] = (*f_up)[i] + (*f_dn)[i];
    }

    return out;
}

} // namespace qss::chem
