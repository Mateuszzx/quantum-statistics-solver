#include <gtest/gtest.h>
#include <cmath>
#include "../../src/chemical_potential/chem_helpers.hpp"
#include "../../include/qss/distributions.hpp"
#include "../../include/qss/types.hpp"

using namespace qss;
using namespace qss::chem;

static Mat make_energy(std::size_t N) {
    Mat E(N, Vec(N));
    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i)
            E[j][i] = (double)i + 0.1 * (double)j;
    return E;
}

TEST(DistributionTotal, DeltaZeroMatchesSingleChannel) {
    const std::size_t N = 8;
    Mat E = make_energy(N);
    double mu = 0.5;
    double kBT = 4.0;
    DistType dist = DistType::FD;

    // Oczekiwany: zwykły rozkład
    auto f_expected = qss::distribution(E[0], mu, kBT, 0.0, 0, 0.0, dist);
    ASSERT_TRUE(f_expected);

    Mat total = distribution_total(E, mu, kBT, dist, 0.0, 0.0);
    // Sprawdź pierwszy wiersz i kilka punktów
    for (std::size_t i = 0; i < N; ++i)
        EXPECT_NEAR(total[0][i], 2.0 * (*f_expected)[i], 1e-12);
}

TEST(DistributionTotal, SpinSplittingEqualsSumOfShiftedChannels) {
    const std::size_t N = 6;
    Mat E = make_energy(N);
    double mu = 0.3;
    double kBT = 3.0;
    double delta = 0.2;
    DistType dist = DistType::FD;

    Mat total = distribution_total(E, mu, kBT, dist, 0.0, delta);

    for (std::size_t j = 0; j < N; ++j) {
        // Przygotuj dwa wektory energii
        Vec up(N), dn(N);
        for (std::size_t i = 0; i < N; ++i) {
            up[i] = E[j][i] - delta;
            dn[i] = E[j][i] + delta;
        }
        auto f_up = qss::distribution(up, mu, kBT, 0.0, 0, 0.0, dist);
        auto f_dn = qss::distribution(dn, mu, kBT, 0.0, 0, 0.0, dist);
        ASSERT_TRUE(f_up && f_dn);
        for (std::size_t i = 0; i < N; ++i) {
            double expected = (*f_up)[i] + (*f_dn)[i];
            EXPECT_NEAR(total[j][i], expected, 1e-12);
        }
    }
}

TEST(DistributionTotal, SymmetryInDelta) {
    const std::size_t N = 5;
    Mat E = make_energy(N);
    double mu = 0.1;
    double kBT = 2.5;
    double delta = 0.35;
    DistType dist = DistType::FD;

    Mat plusD = distribution_total(E, mu, kBT, dist, 0.0, delta);
    Mat minusD = distribution_total(E, mu, kBT, dist, 0.0, -delta);

    for (std::size_t j = 0; j < N; ++j)
        for (std::size_t i = 0; i < N; ++i)
            EXPECT_NEAR(plusD[j][i], minusD[j][i], 1e-12);
}

TEST(DistributionTotal, ThrowsOnRaggedRows) {
    Mat E;
    E.push_back(Vec{0.0, 1.0, 2.0});
    E.push_back(Vec{0.5, 1.5}); // krótszy wiersz

    EXPECT_THROW(
        (void)distribution_total(E, 0.0, 1.0, DistType::FD, 0.0, 0.0),
        std::runtime_error);
}

TEST(DistributionTotal, FiniteValues) {
    const std::size_t N = 7;
    Mat E = make_energy(N);
    double mu = 0.4;
    double kBT = 5.0;
    double delta = 0.15;
    DistType dist = DistType::FD;

    Mat total = distribution_total(E, mu, kBT, dist, 0.0, delta);
    for (auto &row : total)
        for (double v : row)
            EXPECT_TRUE(std::isfinite(v));
}
