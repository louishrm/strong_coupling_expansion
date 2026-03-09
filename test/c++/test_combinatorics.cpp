#include <gtest/gtest.h>
#include "sc_expansion/combinatorics.hpp"
#include <vector>

using namespace sc_expansion;

TEST(CombinatoricsTest, Factorial) {
  EXPECT_EQ(factorial(0), 1);
  EXPECT_EQ(factorial(1), 1);
  EXPECT_EQ(factorial(2), 2);
  EXPECT_EQ(factorial(3), 6);
  EXPECT_EQ(factorial(4), 24);
  EXPECT_EQ(factorial(5), 120);
  EXPECT_EQ(factorial(10), 3628800);
}

TEST(CombinatoricsTest, PartitionGenerator) {
  // Partitions of 4 with max_val 2:
  // {1, 1, 1, 1}
  // {1, 1, 2}
  // {2, 2}
  PartitionGenerator pg(4, 2);
  
  std::vector<std::vector<int>> expected = {
    {1, 1, 1, 1},
    {1, 1, 2},
    {2, 2}
  };
  
  std::vector<std::vector<int>> result;
  while (pg.is_valid()) {
    result.push_back(pg.current());
    if (!pg.next()) break;
  }
  
  EXPECT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < expected.size() && i < result.size(); ++i) {
    EXPECT_EQ(result[i], expected[i]);
  }
}

TEST(CombinatoricsTest, PartitionGeneratorMaxValConstraint) {
  // Partitions of 3 with max_val 1:
  // {1, 1, 1}
  PartitionGenerator pg(3, 1);
  
  EXPECT_TRUE(pg.is_valid());
  EXPECT_EQ(pg.current(), std::vector<int>({1, 1, 1}));
  EXPECT_FALSE(pg.next());
}

TEST(CombinatoricsTest, PartitionGeneratorReset) {
  PartitionGenerator pg(2, 2);
  // {1, 1}, {2}
  
  EXPECT_TRUE(pg.is_valid());
  pg.next();
  pg.next();
  EXPECT_FALSE(pg.is_valid());
  
  pg.reset();
  EXPECT_TRUE(pg.is_valid());
  EXPECT_EQ(pg.current(), std::vector<int>({1, 1}));
}

TEST(CombinatoricsTest, PartitionGeneratorOfSixMaxThree) {
  // Partitions of 6 with max_val 3:
  // {1, 1, 1, 1, 1, 1}
  // {1, 1, 1, 1, 2}
  // {1, 1, 1, 3}
  // {1, 1, 2, 2}
  // {1, 2, 3}
  // {2, 2, 2}
  // {3, 3}
  PartitionGenerator pg(6, 3);
  
  std::vector<std::vector<int>> expected = {
    {1, 1, 1, 1, 1, 1},
    {1, 1, 1, 1, 2},
    {1, 1, 1, 3},
    {1, 1, 2, 2},
    {1, 2, 3},
    {2, 2, 2},
    {3, 3}
  };
  
  std::vector<std::vector<int>> result;
  while (pg.is_valid()) {
    result.push_back(pg.current());
    if (!pg.next()) break;
  }
  
  EXPECT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < expected.size() && i < result.size(); ++i) {
    EXPECT_EQ(result[i], expected[i]);
  }
}

TEST(CombinatoricsTest, PermutationSign) {
  // S4 even permutation: (1, 2, 0, 3) which is (0 1 2)(3)
  EXPECT_DOUBLE_EQ(compute_permutation_sign({1, 2, 0, 3}), 1.0);
  
  // S4 odd permutation: (1, 0, 2, 3) which is (0 1)(2)(3)
  EXPECT_DOUBLE_EQ(compute_permutation_sign({1, 0, 2, 3}), -1.0);
  
  // Identity: (0, 1, 2, 3)
  EXPECT_DOUBLE_EQ(compute_permutation_sign({0, 1, 2, 3}), 1.0);
}

TEST(CombinatoricsTest, SJT_PermutationCount) {
  // n = 1
  SJT sjt1(1);
  int count1 = 1;
  while (sjt1.next_permutation()) { count1++; }
  EXPECT_EQ(count1, 1);

  // n = 3
  SJT sjt3(3);
  int count3 = 1;
  while (sjt3.next_permutation()) { count3++; }
  EXPECT_EQ(count3, 6);

  // n = 4
  SJT sjt4(4);
  int count4 = 1;
  while (sjt4.next_permutation()) { count4++; }
  EXPECT_EQ(count4, 24);
}

TEST(CombinatoricsTest, SJT_PermutationSequence) {
  SJT sjt(3);
  std::vector<std::vector<int>> result;
  result.push_back(sjt.get_permutation());
  while (sjt.next_permutation()) { result.push_back(sjt.get_permutation()); }

  std::vector<std::vector<int>> expected = {{1, 2, 3}, {1, 3, 2}, {3, 1, 2}, {3, 2, 1}, {2, 3, 1}, {2, 1, 3}};

  EXPECT_EQ(result.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i) { EXPECT_EQ(result[i], expected[i]); }
}

