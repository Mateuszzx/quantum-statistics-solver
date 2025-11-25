#include <gtest/gtest.h>
// math utils tests
#include <gtest/gtest.h>
#include <qss/math_utils.hpp>
#include <qss/types.hpp>
#include <cmath>

using namespace qss::math;


TEST(MathUtils, MeanAbs2DZeroForEmpty)
{
    qss::Mat M; // empty
    EXPECT_DOUBLE_EQ(mean_abs_2d(M),0);
}

TEST(MathUtils, MeanAbs2DSimple)
{
    qss::Mat M = {{-1.0, 2.0}, {3.0, -4.0}}; // abs sum = 1+2+3+4 =10
    double N = 2.0;                     // mean_abs_2d uses sum/(N*N)
    double v = mean_abs_2d(M);
    EXPECT_DOUBLE_EQ(v, 10.0 / (N * N));
}
