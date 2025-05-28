//#include "sc_expansion/sc_expansion.hpp"
#include "sc_expansion/hubbard_dimer.hpp"

int main(int argc, char *argv[]) {

  double U    = 1.0;
  double mu   = 0.5;
  double beta = 1.0;

  double n = hubbard_dimer::occupation(U, mu, beta);

  std::cout << "Occupation: " << n << std::endl;

  std::vector<std::pair<int, int>> path = {{0, 2}, {2, 0}};
  std::vector<double> tau               = {0.0, 0.0};
  double G0                             = hubbard_dimer::unperturbed_green_function(U, mu, beta, path, tau);
  std::cout << "Unperturbed Green's function: " << G0 << std::endl;
  return 0;
}
