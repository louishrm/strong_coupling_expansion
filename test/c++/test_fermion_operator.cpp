#include <gtest/gtest.h>
#include "sc_expansion/fermion_operator.hpp"

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

  FermionOperator<1> op = FermionOperator<1>(0); //c_dn
  EXPECT_EQ(op.get_orbital_index(), 0);
  EXPECT_EQ(op.get_action(), 0);

  op = FermionOperator<1>(3); //cdag_up
  EXPECT_EQ(op.get_orbital_index(), 1);
  EXPECT_EQ(op.get_action(), 1);

  FockState<1> state = FockState<1>(1);           //|down>
  auto transition    = op.act_on_state(state);    //cdag_up cdag_down |0> = -|down up>
  EXPECT_EQ(transition.connected_state.state, 3); //|down up>
  EXPECT_DOUBLE_EQ(transition.matrix_element, -1.0);
}
