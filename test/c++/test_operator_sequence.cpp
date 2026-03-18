#include <gtest/gtest.h>
#include "../c++/sc_expansion/operator_sequence.hpp"
#include <vector>

using namespace sc_expansion;

TEST(OperatorSequenceTest, SortingAndSign) {
  // Test time-ordering and fermionic permutation sign.
  // Input: {t1, t2, t3, t4}
  // Let t1 = 0.5, t2 = 0.8, t3 = 0.2, t4 = 0.9
  // ops: {0, 1, 2, 3}
  // Sorted times descending: 0.9, 0.8, 0.5, 0.2 -> corresponding to t4, t2, t1, t3
  // indices: 3, 1, 0, 2
  // Permutation (3, 1, 0, 2) from (0, 1, 2, 3):
  // (0 1 2 3) -> swap 0,3 -> (3 1 2 0) -> swap 1,2 -> (3 2 1 0) -> swap 1,2 -> (3 1 2 0)
  // Let's just trust the compute_permutation_sign, but sign should be deterministic.
  // Inversions for (3, 1, 0, 2): (3,1), (3,0), (3,2), (1,0) -> 4 inversions. Sign = +1.

  std::vector<double> taus = {0.5, 0.8, 0.2, 0.9};
  std::vector<int> ops = {0, 1, 2, 3};

  OperatorSequence seq(taus, ops);

  EXPECT_EQ(seq.taus[0], 0.9);
  EXPECT_EQ(seq.taus[1], 0.8);
  EXPECT_EQ(seq.taus[2], 0.5);
  EXPECT_EQ(seq.taus[3], 0.2);

  EXPECT_EQ(seq.ops[0], 3);
  EXPECT_EQ(seq.ops[1], 1);
  EXPECT_EQ(seq.ops[2], 0);
  EXPECT_EQ(seq.ops[3], 2);

  EXPECT_DOUBLE_EQ(seq.permutation_sign, 1.0);
}

TEST(OperatorSequenceTest, VerifyFlavorAlternation) {
  // Hubbard Atom: Bit 0 is Action, Bit 1 is Spin.
  // We provide already sorted sequences (time is 1.0, 0.9, etc).
  // Flavor 0 (Spin Down): Action 1 (create) = 1, Action 0 (destroy) = 0.
  // Flavor 1 (Spin Up): Action 1 (create) = 3, Action 0 (destroy) = 2.

  // Valid Sequence: cdag_up, c_up, cdag_dn, c_dn (descending time)
  // ops: 3, 2, 1, 0
  OperatorSequence valid_seq({0.4, 0.3, 0.2, 0.1}, {3, 2, 1, 0});
  EXPECT_TRUE(valid_seq.verify_flavor_alternation(0));

  // Invalid: cdag_up, cdag_up, c_up, c_up
  // ops: 3, 3, 2, 2
  OperatorSequence invalid_seq({0.4, 0.3, 0.2, 0.1}, {3, 3, 2, 2});
  EXPECT_FALSE(invalid_seq.verify_flavor_alternation(0));

  // Valid Interleaved: cdag_up, cdag_dn, c_up, c_dn
  // ops: 3, 1, 2, 0
  OperatorSequence interleaved_seq({0.4, 0.3, 0.2, 0.1}, {3, 1, 2, 0});
  EXPECT_TRUE(interleaved_seq.verify_flavor_alternation(0));
}

TEST(OperatorSequenceTest, VerifyInfiniteU) {
  // Infinite U constraint for Hubbard Atom:
  // Global strict alternation of action, and identical spin for (create, destroy) pairs.
  // Since time is descending, first operator in pair is `late` (higher time), second is `early` (lower time).
  // Wait, my verify_infinite_U_constraints rule A:
  // if (early_action == late_action) return false; -> global alternation
  // Rule B: if early is create (1), then late spin == early spin.
  // Let's test this carefully.
  // Late is applied AFTER early. So `early` is applied first (e.g. from empty state).
  // So `early` must be `create`.
  // Wait, time descending: taus[0] > taus[1] > taus[2].
  // Operators are applied from right to left (earliest time first).
  // So ops[n-1] is applied first.

  // Sequence 1: cdag_up(0.4), c_up(0.1)
  // applied: c_up first? No, if empty, we must create first.
  // So ops[1] (time 0.1) is cdag_up (3). ops[0] (time 0.4) is c_up (2).
  // Let's trace it:
  // ops[0] = 2 (destroy up)
  // ops[1] = 3 (create up)
  OperatorSequence valid({0.4, 0.1}, {2, 3});
  EXPECT_TRUE(valid.verify_infinite_U_constraints());

  // Different spin for create and destroy pair
  OperatorSequence invalid1({0.4, 0.1}, {0, 3}); // destroy dn, create up
  EXPECT_FALSE(invalid1.verify_infinite_U_constraints());

  // Same action globally
  OperatorSequence invalid2({0.4, 0.1}, {3, 3});
  EXPECT_FALSE(invalid2.verify_infinite_U_constraints());
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
