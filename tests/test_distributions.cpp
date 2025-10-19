#include <gtest/gtest.h>
#include <qss/distributions.hpp>  // your public API header

TEST(Distributions, FermiDiracBasic) {
  // Replace with your real API. This is a placeholder example:
  const double x = 0.5;
  const double mu = 0.1;
  const double kT = 1.0;
  const double val = qss::fd_distribution(x, mu, kT);
  EXPECT_GE(val, 0.0);
  EXPECT_LE(val, 1.0);
}
