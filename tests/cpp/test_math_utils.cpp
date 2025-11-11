#include <gtest/gtest.h>
// math utils tests
#include <gtest/gtest.h>
#include <qss/math_utils.hpp>
#include <qss/types.hpp>
#include <cmath>

using namespace qss::math;

TEST(MathUtils, MakeEmptyMatDimensions)
{
    auto M = make_empty_mat(3, 5);
    ASSERT_EQ(M.size(), 3u);
    for (const auto &row : M)
        ASSERT_EQ(row.size(), 5u);
}

TEST(MathUtils, MeanAbs2DZeroForEmpty)
{
    qss::Mat M; // empty
    EXPECT_DOUBLE_EQ(mean_abs_2d(M, 10.0), 0.0);
}

TEST(MathUtils, MeanAbs2DSimple)
{
    qss::Mat M = {{-1.0, 2.0}, {3.0, -4.0}}; // abs sum = 1+2+3+4 =10
    double N = 2.0;                     // mean_abs_2d uses sum/(N*N)
    double v = mean_abs_2d(M, N);
    EXPECT_DOUBLE_EQ(v, 10.0 / (N * N));
}

#include <functional>

TEST(MathUtils, FindBracketExpansionSuccess)
{
    auto f = [](double x) { return x - 5.0; };
    double a = 10.0, b = 11.0; // does not bracket root 5
    double outA = 0.0, outB = 0.0;
    bool ok = find_bracket_expansion(f, a, b, 1.5, 25, &outA, &outB);
    ASSERT_TRUE(ok);
    EXPECT_LE(outA, 5.0);
    EXPECT_GE(outB, 5.0);
}

TEST(MathUtils, FindBracketExpansionFail)
{
    auto f = [](double x) { return x * x + 1.0; }; // no real root
    double outA = 0.0, outB = 0.0;
    bool ok = find_bracket_expansion(f, -1.0, 1.0, 1.5, 10, &outA, &outB);
    EXPECT_FALSE(ok);
}
