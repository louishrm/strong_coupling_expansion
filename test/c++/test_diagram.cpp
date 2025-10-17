#include "../../c++/sc_expansion/diagram.hpp"
#include <iostream>

int main() {

  //adjmat adjmat_test = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  adjmat adjmat_test = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}}; //3-cycle with double lines
  Diagram D(adjmat_test);
  bool test1 = D.is_connected();
  bool test2 = D.is_particle_number_conserving();
  int sf     = D.get_symmetry_factor();
  std::cout << "Connected: " << test1 << ", PNC: " << test2 << ", Symmetry factor: " << sf << std::endl;
  return 0;
}