// #include "../../c++/sc_expansion/diagram.hpp"
// #include "../../c++/sc_expansion/hubbard_solver.hpp"
// #include "../../c++/sc_expansion/cumulant.hpp"
// #include <chrono>

// using namespace sc_expansion;

// std::function<std::vector<double>(std::vector<double>)> make_x_to_tau(double beta) {

//   return [beta](std::vector<double> xs) -> std::vector<double> {
//     int n = xs.size();
//     std::vector<double> taus(xs.size());
//     taus[0] = beta * std::pow(xs[0], 1.0 / (double)n);
//     for (int i = 1; i < xs.size(); i++) { taus[i] = taus[i - 1] * std::pow(xs[i], 1.0 / ((double)(n - i))); }
//     return taus;
//   };
// }

// double compute_4th_order_diagrams(std::vector<double> const &taus, std::vector<Diagram> const &diagrams) {

//   double diagram_sum = 0.0;
//   for (auto const &diagram : diagrams) {

//     double val = diagram.evaluate_at_taus(taus, false);
//     diagram_sum += val;
//   }
//   return diagram_sum;
// }

// int main(int argc, char *argv[]) {

//   auto start = std::chrono::high_resolution_clock::now();
//   double x1  = std::stod(argv[1]);
//   double x2  = std::stod(argv[2]);
//   double x3  = std::stod(argv[3]);
//   double x4  = std::stod(argv[4]);

//   double U    = 8.0;
//   double beta = 1.0;
//   double mu   = 2.0;

//   triqs::hilbert_space::fundamental_operator_set fops     = hubbard_atom::make_fops();
//   triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);
//   triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

//   adjmat D4a                    = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
//   adjmat D4b                    = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
//   adjmat D4c                    = {{0, 2}, {2, 0}};                                         //2-cycle with double lines
//   std::vector<Diagram> diagrams = {Diagram(D4a), Diagram(D4b), Diagram(D4c)};

//   for (int i = 0; i < 100; i++) {
//     double diagram_sum = 0.0;
//     for (auto const &diagram : diagrams) {

//       std::vector<double> x = {x1, x2, x3, x4};
//       auto taus             = make_x_to_tau(beta)(x);
//       double val            = diagram.evaluate_at_taus(ad, beta, taus);
//       diagram_sum += val;
//     }
//   }
//   auto end                              = std::chrono::high_resolution_clock::now();
//   std::chrono::duration<double> elapsed = end - start;
//   std::cout << "Elapsed time: " << elapsed.count() << " seconds\n";
//   //std::cout << "Value computed: " << diagram_sum << std::endl;
//   //std::cout.write(reinterpret_cast<const char *>(&diagram_sum), sizeof(diagram_sum));
//   return 0;
// }
