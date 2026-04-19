#include <gtest/gtest.h>
#include <cmath>
#include "../c++/sc_expansion/args.hpp"
#include <vector>

using namespace sc_expansion;

// --- Args Tests ---

TEST(ArgsTest, ConstructorAndSorting) {
  // Scenario: t1 < t2, but input as {t1, t2}. Should be sorted to {t2, t1}.
  // Index 0: spin 0 (down) -> creation (op=2) -- WAIT, N_sites=1, so action bit is bit 1.
  // Action bit is 1 << N_sites. For N_sites=1, action bit is 2.
  // op: bit 0 is orbital (0=down, 1=up), bit 1 is action (0=destroy, 1=create).
  // annihilate down: 0
  // annihilate up:   1
  // create down:     2
  // create up:       3

  std::vector<double> taus = {0.2, 0.8};
  std::vector<FermionOperator<double>> ops = {
    FermionOperator<double>(3), // create up
    FermionOperator<double>(1)  // annihilate up
  };

  Args<double> args(taus, ops);

  EXPECT_EQ(args.order, 2);
  // Sorted descending:
  EXPECT_DOUBLE_EQ(args.taus[0], 0.8);
  EXPECT_DOUBLE_EQ(args.taus[1], 0.2);

  // After sorting (argsort={1, 0}):
  EXPECT_EQ(args.ops[0].op, 1); // annihilate up (from tau=0.8)
  EXPECT_EQ(args.ops[1].op, 3); // create up (from tau=0.2)

  // One swap: permutation sign should be -1.0
  EXPECT_DOUBLE_EQ(args.permutation_sign, -1.0);
}

TEST(ArgsTest, PermutationSignComplex) {
  // 4 operators, multiple swaps
  std::vector<double> taus = {0.1, 0.4, 0.2, 0.3}; // argsort: {1, 3, 2, 0}
  std::vector<FermionOperator<double>> ops = {
    FermionOperator<double>(0),
    FermionOperator<double>(2),
    FermionOperator<double>(1),
    FermionOperator<double>(3)
  };
  Args<double> args(taus, ops);

  // inversions in {1, 3, 2, 0}:
  // (1,0), (3,2), (3,0), (2,0) -> 4 inversions. Sign = +1.0.
  EXPECT_DOUBLE_EQ(args.permutation_sign, 1.0);
}

TEST(ArgsTest, VerifyConsecutiveTermsInfiniteU) {
  // Test valid sequence: cup(0.8), cdag_up(0.2)
  // sorted: cup(0.8), cdag_up(0.2)
  // ops: 1, 3
  // sequence: 1(annihilate up), 3(create up)
  // Path: |up> --(1)--> |0> --(3)--> |up>  VALID
  {
    std::vector<double> taus = {0.8, 0.2};
    std::vector<FermionOperator<double>> ops = {
        FermionOperator<double>(1), // cup
        FermionOperator<double>(3)  // cdag_up
    };
    Args<double> args(taus, ops);
    EXPECT_TRUE(args.operator_sequence_is_valid_infinite_U());
  }

  // Invalid: Consecutive creations on same orbital (Pauli exclusion)
  {
    std::vector<double> taus = {0.8, 0.2};
    std::vector<FermionOperator<double>> ops = {
        FermionOperator<double>(3), // cdag_up
        FermionOperator<double>(3)  // cdag_up
    };
    Args<double> args(taus, ops);
    EXPECT_FALSE(args.operator_sequence_is_valid_infinite_U());
  }

  // Invalid: Spin mismatch (creation of up, destruction of down)
  // sorted: cdag_up(0.8), cdn(0.2)
  // sequence: 3(create up), 0(annihilate down)
  // start |0>: |0> --(3)--> |up> --(0)--> INVALID (no down to annihilate)
  // start |down>: |down> --(3)--> |up down> (forbidden)
  {
    std::vector<double> taus = {0.8, 0.2};
    std::vector<FermionOperator<double>> ops = {
        FermionOperator<double>(3), // cdag_up
        FermionOperator<double>(0)  // cdn
    };
    Args<double> args(taus, ops);
    EXPECT_FALSE(args.operator_sequence_is_valid_infinite_U());
  }

  // Valid: cup(0.9), cdag_up(0.7), cdn(0.5), cdag_dn(0.3)
  // sorted: cup(0.9), cdag_up(0.7), cdn(0.5), cdag_dn(0.3)
  // sequence: 1, 3, 0, 2
  // start |up>: |up> --(1)--> |0> --(3)--> |up> --(0)--> INVALID (no down)
  // wait... |up> --(1)--> |0> --(3)--> |up>  ... then 0, 2?
  // Let's try start |up down>? NO, forbidden.
  // How about: cup(0.9), cdag_up(0.7), cdn(0.5), cdag_dn(0.3)?
  // If we start at |up>, we annihilate up at 0.9, then create up at 0.7.
  // Then we need to handle down. If we had down, we could annihilate it at 0.5 and create it at 0.3.
  // So start state |up down>? NO.
  // The only valid start states are |0>, |up>, |down>.
  // If we start at |up>, we can do 1, 3. But then we are back at |up>, and we need to do 0, 2. Still no down.
  // If we start at |0>, we can't do 1.
  
  // Let's try: cup(0.9), cdag_up(0.7), cdn(0.5), cdag_dn(0.3)
  // This requires starting with both up and down? NO.
  // Wait, the operators are in a trace.
  // Tr(rho A B C D) = Tr(rho c_up(0.9) cdag_up(0.7) c_dn(0.5) cdag_dn(0.3))
  // This is only non-zero if rho has components that allow this.
  // But infinite U forbids |up down>.
  // So rho can only be |0><0|, |up><up|, |down><down|.
  // None of these allow both up and down operators in a single term of the trace
  // UNLESS they don't overlap in time such that |up down| is never reached.
  // 0.9: cup -> needs |up> or |up down|.
  // 0.7: cdag_up -> needs |0> or |down|.
  // 0.5: cdn -> needs |down> or |up down|.
  // 0.3: cdag_dn -> needs |0> or |up|.
  
  // Trace path:
  // Start |up>: |up> --(0.9: cup)--> |0> --(0.7: cdag_up)--> |up> --(0.5: cdn)--> INVALID
  // Start |down>: |down> --(0.9: cup)--> INVALID
  // Start |0>: |0> --(0.9: cup)--> INVALID
  
  // So the sequence 1, 3, 0, 2 is actually INVALID for infinite U!
  // A valid one with 4 ops would be: 1, 3, 1, 3 or something.
  // Or 1, 3, then later 0, 2? No, they are all in one trace.
  // 1, 3, 0, 2 is valid for FINITE U.
  
  // Let's use a truly valid one for infinite U: cup(0.9), cdag_up(0.7), cup(0.5), cdag_up(0.3)
  // Start |up>: |up> --(1)--> |0> --(3)--> |up> --(1)--> |0> --(3)--> |up>  VALID!
  {
    std::vector<double> taus = {0.9, 0.7, 0.5, 0.3};
    std::vector<FermionOperator<double>> ops = {
        FermionOperator<double>(1), // cup
        FermionOperator<double>(3), // cdag_up
        FermionOperator<double>(1), // cup
        FermionOperator<double>(3)  // cdag_up
    };
    Args<double> args(taus, ops);
    EXPECT_TRUE(args.operator_sequence_is_valid_infinite_U());
  }
}
