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
