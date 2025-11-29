#pragma once
#include <stdexcept>
#include "qss/distributions.hpp"
#include "qss/types.hpp"


namespace qss::chem {

    constexpr double SPIN_UP_SHIFT = 1.0;
    constexpr double SPIN_DN_SHIFT = -1.0;

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

        auto spin_up = qss::distribution(E[j], mu, beta, delta, static_cast<int>(SPIN_UP_SHIFT), dist);
        auto spin_dn = qss::distribution(E[j], mu, beta, delta, static_cast<int>(SPIN_DN_SHIFT), dist);
        if (!spin_up || !spin_dn)
            throw std::runtime_error("distribution() failed");

        // No splitting: copy f(E) results to output
        for (std::size_t i = 0; i < nx; ++i) {
            out[j][i] = (*spin_up)[i] + (*spin_dn)[i];
        }
    }

    return out;
}

} // namespace qss::chem
