// src/chem_helpers.hpp
#pragma once
#include "qss/distributions.hpp"
#include "qss/types.hpp"
#include <stdexcept> // add this to fix std::runtime_error errors

namespace qss::chem
{

// Builds total distribution map for one mu/beta:
//   if delta==0 -> F = f(E)
//   else        -> F = f(E - delta) + f(E + delta)  (↑/↓ channels)
// Returns Mat with shape E.shape()
inline Mat distribution_total(const Mat &E,
                              double mu,
                              double beta,
                              DistType dist,
                              double delta)
{
    const std::size_t ny = E.size();
    const std::size_t nx = (ny ? E[0].size() : 0);
    Mat out(ny, Vec(nx));

    for (std::size_t j = 0; j < ny; ++j)
    {
        // base channel
        auto f0 = qss::distribution(E[j], mu, beta, dist);
        if (!f0)
            throw std::runtime_error("distribution() failed");
        if (delta == 0.0)
        {
            out[j] = std::move(*f0);
            continue;
        }
        // spin split: up/down energies
        Vec row_up(nx), row_dn(nx);
        for (std::size_t i = 0; i < nx; ++i)
        {
            row_up[i] = E[j][i] - delta;
            row_dn[i] = E[j][i] + delta;
        }
        auto f_up = qss::distribution(row_up, mu, beta, dist);
        auto f_dn = qss::distribution(row_dn, mu, beta, dist);
        if (!f_up || !f_dn)
            throw std::runtime_error("distribution() failed (split)");
        for (std::size_t i = 0; i < nx; ++i)
            out[j][i] = (*f_up)[i] + (*f_dn)[i];
    }
    return out;
}

} // namespace qss::chem
