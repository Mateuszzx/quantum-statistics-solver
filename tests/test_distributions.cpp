// tests/test_distributions.cpp
#include <gtest/gtest.h>
#include <qss/distributions.hpp>
#include <qss/types.hpp>          // for UP_EXP_LIMIT / DOWN_EXP_LIMIT (public)

using qss::DistType;

TEST(Distributions, InvalidEnumReturnsNullopt) {
  // Pretend-cast an invalid enum (not in DistType)
  const auto bad = static_cast<DistType>(999);
  auto o1 = qss::distribution(0.0, 0.0, 1.0, bad);
  auto o2 = qss::distribution(qss::Vec{0.0, 1.0}, 0.0, 1.0, bad);
  EXPECT_FALSE(o1.has_value());
  EXPECT_FALSE(o2.has_value());
}

TEST(Distributions, VectorEmptyOK) {
  const qss::Vec empty{};
  auto vopt = qss::distribution(empty, 0.0, 1.0, DistType::FD);
  ASSERT_TRUE(vopt.has_value());
  EXPECT_TRUE(vopt->empty());
}



TEST(Distributions, FD_Scalar_Range) {
  const double mu = 0.1;
  const double beta = 1.0;                 // kT = 1 -> beta = 1
  const double ek  = 0.5;

  auto opt = qss::distribution(ek, mu, beta, DistType::FD);
  ASSERT_TRUE(opt.has_value());
  const double v = *opt;

  // With spin degeneracy 2.0, FD in [0, 2]
  EXPECT_GE(v, 0.0);
  EXPECT_LE(v, 2.0);
}

TEST(Distributions, FD_Scalar_Limits_Saturate) {
  const double mu = 0.0;
  const double beta = 1.0;

  // x = beta*(ek - mu)
  // Large positive x -> ~0
  {
    const double ek_pos = qss::UP_EXP_LIMIT + 100.0; // definitely > limit
    auto opt = qss::distribution(ek_pos, mu, beta, DistType::FD);
    ASSERT_TRUE(opt.has_value());
    EXPECT_NEAR(*opt, 0.0, 0.0);
  }

  // Large negative x -> ~2 (spin=2)
  {
    const double ek_neg = qss::DOWN_EXP_LIMIT - 100.0; // definitely < limit
    auto opt = qss::distribution(ek_neg, mu, beta, DistType::FD);
    ASSERT_TRUE(opt.has_value());
    EXPECT_NEAR(*opt, 2.0, 0.0);
  }
}

TEST(Distributions, FD_Scalar_MonotonicInEnergy) {
  const double mu = 0.0;
  const double beta = 2.0;  // steeper curve

  const double e_low  = -1.0;
  const double e_high =  +1.0;

  const double v_low  = *qss::distribution(e_low,  mu, beta, DistType::FD);
  const double v_high = *qss::distribution(e_high, mu, beta, DistType::FD);

  // As energy increases (with fixed mu, beta>0), FD occupancy decreases
  EXPECT_GT(v_low, v_high);
}

TEST(Distributions, FD_BetaZero_LimitsToHalfSpin) {
  const double mu = 0.0;
  const double beta = 0.0;   // beta -> 0 => x = 0, logistic = 1/2, spin=2 -> value = 1.0
  const double ek = 123.456;

  const double v = *qss::distribution(ek, mu, beta, DistType::FD);
  EXPECT_NEAR(v, 1.0, 1e-12);  // half-occupancy with spin=2
}

TEST(Distributions, FD_Vector_Equals_ElementwiseScalar) {
  const double mu = 0.2;
  const double beta = 1.5;
  const qss::Vec ek = {-5.0, -1.0, 0.0, 0.3, 2.5, 10.0};

  auto vopt = qss::distribution(ek, mu, beta, DistType::FD);
  ASSERT_TRUE(vopt.has_value());
  const qss::Vec& v = *vopt;

  ASSERT_EQ(v.size(), ek.size());
  for (std::size_t i = 0; i < ek.size(); ++i) {
    const double s = *qss::distribution(ek[i], mu, beta, DistType::FD);
    EXPECT_NEAR(v[i], s, 1e-12) << "Mismatch at i=" << i;
  }
}
