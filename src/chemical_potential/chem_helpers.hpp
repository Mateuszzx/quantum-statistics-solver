// src/chem_helpers.hpp
#pragma once
#include <stdexcept>
#include "qss/distributions.hpp"
#include "qss/types.hpp"

namespace qss::chem {

/*
  distribution_total:
  - Buduje całkowitą obsadę (sumę po spinach) dla zadanych:
    E    : siatka energii 2D (ny wierszy x nx kolumn)
    mu   : potencjał chemiczny
    beta : odwrotność temperatury
    dist : typ rozkładu (FD/FLB/NFLB/...)
    delta: przesunięcie Zeemana (g*muB*H/2) w jednostkach E

  Zasada:
  - Dla delta == 0: brak rozszczepienia spinowego → F = f(E)
  - Dla delta != 0: rozszczepienie spinowe → F = f(E - delta) + f(E + delta)
    (sumujemy kanały ↑ i ↓ punkt-po-punkcie)
*/
inline Mat distribution_total(const Mat& E, double mu, double beta,
                            DistType dist, double delta) {
    // Wymiary macierzy energii
    const std::size_t ny = E.size();
    const std::size_t nx = ny ? E[0].size() : 0;

    // Prosty test spójności wierszy (opcjonalnie; O(ny))
    for (std::size_t j = 1; j < ny; ++j) {
        if (E[j].size() != nx)
            throw std::runtime_error("distribution_total: E has ragged rows");
    }

    // Wyjście ma taki sam kształt jak E
    Mat out(ny, Vec(nx));

    for (std::size_t j = 0; j < ny; ++j) {
        // Najpierw licz rozkład dla wiersza E[j] bez przesunięcia
        // (przydaje się w przypadku delta == 0)
        auto f0 = qss::distribution(E[j], mu, beta, dist);
        if (!f0)
            throw std::runtime_error("distribution() failed");

        if (delta == 0.0) {
            // Brak rozszczepienia: kopiujemy wyniki f(E) do wyjścia
            out[j] = std::move(*f0);
            continue; // kolejny wiersz
        }

        // Rozszczepienie spinowe: przygotuj dwa „wiersze energii”
        // przesunięte o ±delta i oblicz rozkłady dla obu kanałów.
        Vec row_up(nx), row_dn(nx);
        for (std::size_t i = 0; i < nx; ++i) {
            row_up[i] = E[j][i] - delta; // kanał ↑
            row_dn[i] = E[j][i] + delta; // kanał ↓
        }

        auto f_up = qss::distribution(row_up, mu, beta, dist);
        auto f_dn = qss::distribution(row_dn, mu, beta, dist);
        if (!f_up || !f_dn)
            throw std::runtime_error("distribution() failed (split)");

        // Suma kanałów spinowych punkt-po-punkcie
        for (std::size_t i = 0; i < nx; ++i)
            out[j][i] = (*f_up)[i] + (*f_dn)[i];
    }

    return out;
}

} // namespace qss::chem
