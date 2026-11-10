#pragma once
#include <stdexcept>
#include <algorithm>
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
 * @param kBT Temperature in energy units.
 * @param dist Distribution type (FD/FLB/NFLB/...).
 * @param U double filling parameter
 * @param delta Zeeman shift (g*muB*H/2) in energy units.
 * @return Mat Total distribution matrix.
 *
 * Principle:
 * - For delta == 0: No spin splitting -> F = f(E)
 * - For delta != 0: Spin splitting -> F = f(E - delta) + f(E + delta)
 *   (summing up and down channels point-by-point)
 */
inline Mat distribution_total(const Mat& E, double mu, double kBT,
                            DistType dist, double U, double delta) {
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
    double result;
    for (std::size_t j = 0; j < ny; ++j) {

        auto spin_up = qss::distribution(E[j], mu, kBT, delta, static_cast<int>(SPIN_UP_SHIFT), U, dist);
        auto spin_dn = qss::distribution(E[j], mu, kBT, delta, static_cast<int>(SPIN_DN_SHIFT), U, dist);
        // auto n_kd = std::optional<Vec>{}; // Left empty by the inner redeclaration.
        if (!spin_up || !spin_dn)
            throw std::runtime_error("distribution() failed");
        // if (dist == DistType::HUBBARD){
            // auto n_kd = qss::distribution(E[j], mu, kBT, delta, static_cast<int>(SPIN_UP_SHIFT), U, dist, true); // Shadows the outer n_kd with a new local variable.
        // }
        // No splitting: copy f(E) results to output
        for (std::size_t i = 0; i < nx; ++i) {
            result = std::clamp((*spin_up)[i], 0.0, 1.0) + std::clamp((*spin_dn)[i], 0.0, 1.0);
            if (dist == DistType::FD  || dist == DistType::FLB || dist == DistType::HUBBARD)
            {
                out[j][i] = result;
            }
            // else if (dist == DistType::HUBBARD) // Unreachable: HUBBARD is already handled above.
            // {
                // out[j][i] = result - std::clamp((*n_kd)[i], 0.0, 1.0); // Subtracting doublons undercounts electrons.
            // }
            else
            {
                out[j][i] = std::clamp(result, 0.0, 1.0);
            }
        }
    }

    return out;
}

} // namespace qss::chem
