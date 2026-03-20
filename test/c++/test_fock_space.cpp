#include <gtest/gtest.h>
#include "sc_expansion/fock_space.hpp"

using namespace sc_expansion;

//--- FockState Tests ---
TEST(FockStateTest, AtomOccupation) {

  FockState<1> state = FockState<1>(0); //empty
  EXPECT_FALSE(state.is_occupied(0));
  EXPECT_FALSE(state.is_occupied(1));

  state = FockState<1>(1); //down
  EXPECT_TRUE(state.is_occupied(0));
  EXPECT_FALSE(state.is_occupied(1));

  state = FockState<1>(2); //up
  EXPECT_FALSE(state.is_occupied(0));
  EXPECT_TRUE(state.is_occupied(1));

  state = FockState<1>(3); //down up
  EXPECT_TRUE(state.is_occupied(0));
  EXPECT_TRUE(state.is_occupied(1));
}

TEST(FockStateTest, DimerOccupation) {

  FockState<2> state = FockState<2>(0); //empty
  EXPECT_FALSE(state.is_occupied(0));
  EXPECT_FALSE(state.is_occupied(1));
  EXPECT_FALSE(state.is_occupied(2));
  EXPECT_FALSE(state.is_occupied(3));

  state = FockState<2>(5); //|down up, 0>
  EXPECT_TRUE(state.is_occupied(0));
  EXPECT_FALSE(state.is_occupied(1));
  EXPECT_TRUE(state.is_occupied(2));
  EXPECT_FALSE(state.is_occupied(3));
}

//--- FermionOperator Tests ---

TEST(FermionOperatorTest, AtomicStates) {

  FermionOperator<1, double> op = FermionOperator<1, double>(0); //c_dn
  EXPECT_EQ(op.get_orbital_index(), 0);
  EXPECT_EQ(op.get_action(), 0);

  op = FermionOperator<1, double>(3); //cdag_up
  EXPECT_EQ(op.get_orbital_index(), 1);
  EXPECT_EQ(op.get_action(), 1);

  FockState<1> state = FockState<1>(1);        //|down>
  auto transition    = op.act_on_state(state); //cdag_up cdag_down |0> = -|down up>
  EXPECT_EQ(transition.connected_state, 3);    //|down up>
  EXPECT_DOUBLE_EQ(transition.matrix_element, -1.0);
}

TEST(FermionOperatorTest, DimerStates) {

  // Basis order: [up_N-1 ... up_0 | down_N-1 ... down_0]
  // Indices for N=2: [3 (2u), 2 (1u) | 1 (2d), 0 (1d)]

  FermionOperator<2, double> cdag_2u = FermionOperator<2, double>((1 << 2) | 3);
  FermionOperator<2, double> c_1u    = FermionOperator<2, double>((0 << 2) | 2);

  // 1. cdag_2u on |1u> (state 4, bit 2)
  FockState<2> state_1u = FockState<2>(4);
  auto trans1           = cdag_2u.act_on_state(state_1u);
  // Expect |2u, 1u> (state 12, bits 2 and 3 set)
  // Sign: 2u jumps over 1u (bit 2), so -1.0
  EXPECT_EQ(trans1.connected_state, 12);
  EXPECT_DOUBLE_EQ(trans1.matrix_element, -1.0);

  // 2. c_1u on |1u, 1d> (state 5, bits 2 and 0 set)
  FockState<2> state_1u1d = FockState<2>(5);
  auto trans2             = c_1u.act_on_state(state_1u1d);
  // Expect |1d> (state 1, bit 0 set)
  // Sign: 1u jumps over 1d (bit 0), so -1.0
  EXPECT_EQ(trans2.connected_state, 1);
  EXPECT_DOUBLE_EQ(trans2.matrix_element, -1.0);

  // 3. Pauli exclusion: cdag_2u on |2u, 1u>
  FockState<2> state_2u1u = FockState<2>(12);
  auto trans3             = cdag_2u.act_on_state(state_2u1u);
  EXPECT_EQ(trans3.connected_state, -1);
  EXPECT_DOUBLE_EQ(trans3.matrix_element, 0.0);
}

TEST(FermionOperatorTest, SparseMatrixAtomic) {
  FermionOperator<1, double> cdag_up = FermionOperator<1, double>((1 << 1) | 1);

  // Create simple eigenstates for the atomic case
  // |0>, |down>, |up>, |up down>
  std::vector<Eigenstate<double>> eigenstates(4);
  eigenstates[0].coefficients = {{0, 1.0}};
  eigenstates[1].coefficients = {{1, 1.0}};
  eigenstates[2].coefficients = {{2, 1.0}};
  eigenstates[3].coefficients = {{3, 1.0}};

  auto matrix = cdag_up.compute_sparse_matrix(eigenstates);

  // cdag_up |0> = |up>  => entry (2, 0) value 1.0
  // cdag_up |down> = -|down up> => entry (3, 1) value -1.0
  // Basis: 0:|0>, 1:|down>, 2:|up>, 3:|up down>

  EXPECT_EQ(matrix.entries.size(), 2);

  bool found1 = false;
  bool found2 = false;
  for (auto const &entry : matrix.entries) {
    if (entry.row == 2 && entry.col == 0) {
      EXPECT_DOUBLE_EQ(entry.value, 1.0);
      found1 = true;
    }
    if (entry.row == 3 && entry.col == 1) {
      EXPECT_DOUBLE_EQ(entry.value, -1.0);
      found2 = true;
    }
  }
  EXPECT_TRUE(found1);
  EXPECT_TRUE(found2);
}
