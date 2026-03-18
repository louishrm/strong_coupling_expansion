#include <gtest/gtest.h>
#include "sc_expansion/dual.hpp"

TEST(DualTest, Constructor) {
  Dual d1;
  EXPECT_DOUBLE_EQ(d1.value, 0.0);
  EXPECT_DOUBLE_EQ(d1.derivative, 0.0);

  Dual d2(1.5, 2.0);
  EXPECT_DOUBLE_EQ(d2.value, 1.5);
  EXPECT_DOUBLE_EQ(d2.derivative, 2.0);
}

TEST(DualTest, Addition) {
  Dual d1(1.0, 0.5);
  Dual d2(2.0, 1.5);
  Dual res = d1 + d2;
  EXPECT_DOUBLE_EQ(res.value, 3.0);
  EXPECT_DOUBLE_EQ(res.derivative, 2.0);
}

TEST(DualTest, Subtraction) {
  Dual d1(5.0, 3.0);
  Dual d2(2.0, 1.0);
  Dual res = d1 - d2;
  EXPECT_DOUBLE_EQ(res.value, 3.0);
  EXPECT_DOUBLE_EQ(res.derivative, 2.0);
}

TEST(DualTest, Multiplication) {
  // (u*v)' = u*v' + u'*v
  Dual d1(2.0, 0.5);
  Dual d2(3.0, 1.0);
  Dual res = d1 * d2;
  EXPECT_DOUBLE_EQ(res.value, 6.0);
  // 2.0 * 1.0 + 0.5 * 3.0 = 2.0 + 1.5 = 3.5
  EXPECT_DOUBLE_EQ(res.derivative, 3.5);
}

TEST(DualTest, Division) {
  // (u/v)' = (u'v - uv') / v^2
  Dual d1(6.0, 2.0);
  Dual d2(2.0, 1.0);
  Dual res = d1 / d2;
  EXPECT_DOUBLE_EQ(res.value, 3.0);
  // (2.0 * 2.0 - 6.0 * 1.0) / 2.0^2 = (4.0 - 6.0) / 4.0 = -0.5
  EXPECT_DOUBLE_EQ(res.derivative, -0.5);
}

TEST(DualTest, UnaryMinus) {
  Dual d1(1.5, -2.0);
  Dual res = -d1;
  EXPECT_DOUBLE_EQ(res.value, -1.5);
  EXPECT_DOUBLE_EQ(res.derivative, 2.0);
}

TEST(DualTest, Exp) {
  // (exp(u))' = exp(u) * u'
  Dual d1(1.0, 2.0);
  Dual res = exp(d1);
  double expected_val = std::exp(1.0);
  EXPECT_DOUBLE_EQ(res.value, expected_val);
  EXPECT_DOUBLE_EQ(res.derivative, expected_val * 2.0);
}

TEST(DualTest, Sqrt) {
  // (sqrt(u))' = 0.5 * u' / sqrt(u)
  Dual d1(4.0, 2.0);
  Dual res = sqrt(d1);
  EXPECT_DOUBLE_EQ(res.value, 2.0);
  // 0.5 * 2.0 / 2.0 = 0.5
  EXPECT_DOUBLE_EQ(res.derivative, 0.5);
}

TEST(DualTest, Power) {
  // (u^p)' = p * u^(p-1) * u'
  Dual d1(2.0, 3.0);
  double p = 3.0;
  Dual res = pow(d1, p);
  EXPECT_DOUBLE_EQ(res.value, 8.0);
  // 3.0 * 2.0^2 * 3.0 = 3.0 * 4.0 * 3.0 = 36.0
  EXPECT_DOUBLE_EQ(res.derivative, 36.0);
}

TEST(DualTest, ScalarOperations) {
  Dual d1(2.0, 0.5);
  double scalar = 3.0;

  // Dual * scalar
  Dual res1 = d1 * scalar;
  EXPECT_DOUBLE_EQ(res1.value, 6.0);
  EXPECT_DOUBLE_EQ(res1.derivative, 1.5);

  // scalar * Dual
  Dual res2 = scalar * d1;
  EXPECT_DOUBLE_EQ(res2.value, 6.0);
  EXPECT_DOUBLE_EQ(res2.derivative, 1.5);

  // Dual / scalar
  Dual res3 = d1 / scalar;
  EXPECT_DOUBLE_EQ(res3.value, 2.0 / 3.0);
  EXPECT_DOUBLE_EQ(res3.derivative, 0.5 / 3.0);
}
