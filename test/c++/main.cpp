//#include "sc_expansion/sc_expansion.hpp"
#include "../c++/sc_expansion/hubbard_atom.hpp"

#include <cmath>
#include <iostream>
#include <xfac/tensor/tensor_ci_2.h>
#include <nda/nda.hpp>
#include "../c++/sc_expansion/cumulant.hpp"
// #include <function>

double dimer_Omega2a(auto ad, double beta, std::vector<double> tau) {

  double G0_12 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta); //G(1|2)
  double G0_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)

  double sign              = 1.0;
  double symmetry_factor   = 1 / 2.0;
  double free_multiplicity = 2.0;
  double spin_factor       = 2.0; // for spin degeneracy
  double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

  return prefactor * G0_12 * G0_21;
}

double dimer_Omega2(auto ad, double beta, double delta) {

  double result = 0.0;

  for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {

    std::vector<double> tau = {0, tau2};

    double result_a = dimer_Omega2a(ad, beta, tau);
    result += result_a * delta;
  }
  return result;
}

double dimer_Omega4a(auto ad, double beta, std::vector<double> tau) {

  double G01_14 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[3], 0}}, ad, beta); //G(1|4)
  double G01_21 = compute_cumulant_decomposition({{tau[1], 0}}, {{tau[0], 0}}, ad, beta); //G(2|1)
  double G01_32 = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[1], 0}}, ad, beta); //G(3|2)
  double G01_43 = compute_cumulant_decomposition({{tau[3], 0}}, {{tau[2], 0}}, ad, beta); //G(4|3)

  double sign              = -1.0;
  double symmetry_factor   = 1 / 4.0;
  double free_multiplicity = 2.0;
  double spin_factor       = 2.0; // for spin degeneracy
  double prefactor         = sign * symmetry_factor * free_multiplicity * spin_factor;

  return prefactor * G01_14 * G01_21 * G01_32 * G01_43;
}

double dimer_Omega4b(auto ad, double beta, std::vector<double> tau) {
  double result = 0.0;

  double C02_up = compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 0}}, {{tau[0], 0}, {tau[2], 0}}, ad, beta); //C(2 4 | 1 3)

  double G01_14_up = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta);
  double G01_32_up = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[3], 0}}, ad, beta);
  result += 2 * C02_up * G01_14_up * G01_32_up;

  double C02_updown = compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 1}}, {{tau[0], 0}, {tau[2], 1}}, ad, beta); //C(2 4 | 1 3)

  double G01_14_up_2 = compute_cumulant_decomposition({{tau[0], 0}}, {{tau[1], 0}}, ad, beta); //G(1|4)
  double G01_32_down = compute_cumulant_decomposition({{tau[2], 1}}, {{tau[3], 1}}, ad, beta); //G(3|2)

  result += 2 * C02_updown * G01_14_up_2 * G01_32_down;

  double sign              = 1.0;
  double symmetry_factor   = 1.0 / 2.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega4c(auto ad, double beta, std::vector<double> tau) {

  double result = 0.0;

  result += 2 * compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 0}}, {{tau[0], 0}, {tau[2], 0}}, ad, beta)
     * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 0}}, {{tau[1], 0}, {tau[3], 0}}, ad, beta);

  result += 2 * compute_cumulant_decomposition({{tau[1], 0}, {tau[3], 1}}, {{tau[0], 0}, {tau[2], 1}}, ad, beta)
     * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 1}}, {{tau[1], 0}, {tau[3], 1}}, ad, beta);

  result += 2 * compute_cumulant_decomposition({{tau[1], 1}, {tau[3], 0}}, {{tau[0], 0}, {tau[2], 1}}, ad, beta)
     * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 1}}, {{tau[1], 1}, {tau[3], 0}}, ad, beta);

  double sign              = 1.0;
  double symmetry_factor   = 1 / 8.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega4(auto ad, double beta, std::vector<double> tau) {
  return dimer_Omega4a(ad, beta, tau) + dimer_Omega4b(ad, beta, tau) + dimer_Omega4c(ad, beta, tau);
}

std::vector<double> linspace(double start, double end, int num, bool endpoint) {
  std::vector<double> result;
  if (num <= 0) { return result; }
  if (num == 1) {
    result.push_back(start);
    return result;
  }
  double step = (end - start) / (endpoint ? (num - 1) : num);
  for (int i = 0; i < num; i++) { result.push_back(start + i * step); }
  return result;
}

std::function<std::vector<double>(std::vector<double>)> make_x_to_tau(double beta) {

  return [beta](std::vector<double> xs) -> std::vector<double> {
    int n = xs.size();
    std::vector<double> taus(xs.size());
    taus[0] = beta * std::pow(xs[0], 1.0 / (double)n);
    for (int i = 1; i < xs.size(); i++) { taus[i] = taus[i - 1] * std::pow(xs[i], 1.0 / ((double)(n - i))); }
    return taus;
  };
}

double sum_TCI(xfac::CTensorCI2<double, double> &TT, std::vector<std::vector<double>> &grid) {

  int n = TT.len();
  arma::mat phi_prev(1, 1, arma::fill::ones);

  for (int i = 0; i < n - 1; i++) {
    arma::mat phi(1, TT.getPivotsAt(i).size(), arma::fill::zeros);
    for (auto index : grid[i]) {
      auto TP1 = TT.get_TP1_at(i, index);
      // std::cout << "TP1 size " << TP1.n_rows << " x " << TP1.n_cols << std::endl;
      // std::cout << "phi_prev size " << phi_prev.n_rows << " x " << phi_prev.n_cols << std::endl;
      phi += phi_prev * TP1;
    }
    phi_prev = phi;
  }

  arma::mat phi(1, 1, arma::fill::zeros);
  for (auto index : grid[n - 1]) {
    auto T = TT.get_T_at(n - 1, index);
    phi += phi_prev * T;
  }
  return phi(0, 0);
}

int main(int argc, char *argv[]) {

  int grid_size = std::stoi(argv[1]);
  int order     = std::stoi(argv[2]);
  int bond_dim  = std::stoi(argv[3]);

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops = hubbard_atom::make_fops();

  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);

  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  std::vector<std::vector<double>> grid;
  for (int i = 0; i < order; i++) {
    auto row = linspace(0.0, 1.0, grid_size + 2, true);
    row.pop_back();         // remove the endpoint 1.0
    row.erase(row.begin()); // remove the start point 0.0
    grid.push_back(row);
  }

  std::vector<double> weights;
  for (int i = 0; i < grid_size; i++) weights.push_back(1.0 / (grid_size));

  // double sum = 0.0;
  // for (auto [wx, wy] : itertools::product(itertools::zip(grid[0], weights), itertools::zip(grid[1], weights))) {
  //   auto [x, w_x] = wx;
  //   auto [y, w_y] = wy;
  //   auto taus     = make_x_to_tau(beta)({x, y});
  //   sum += dimer_Omega2a(ad, beta, taus) * w_x * w_y;
  // }

  double sum = 0.0;
  for (auto [wx1, wx2, wx3, wx4] : itertools::product(itertools::zip(grid[0], weights), itertools::zip(grid[1], weights),
                                                      itertools::zip(grid[2], weights), itertools::zip(grid[3], weights))) {
    auto [x1, w_x1] = wx1;
    auto [x2, w_x2] = wx2;
    auto [x3, w_x3] = wx3;
    auto [x4, w_x4] = wx4;

    auto taus = make_x_to_tau(beta)({x1, x2, x3, x4});
    sum += dimer_Omega4(ad, beta, taus) * w_x1 * w_x2 * w_x3 * w_x4;
  }

  std::cout << "Result from direct summation: " << sum * std::pow(beta, order) << std::endl;

  auto omega2a_adaptator = [&beta, &ad, &weights, &grid](std::vector<double> xs) -> double {
    std::vector<double> taus = make_x_to_tau(beta)(xs);
    double weight            = 1.0;
    for (auto [i, x_i] : itertools::enumerate(xs)) {
      for (auto [j, x_j] : itertools::enumerate(grid[i])) {
        if (x_i == x_j) {
          weight *= weights[j];
          break;
        }
      }
    }
    return dimer_Omega2a(ad, beta, taus) * weight;
  };

  auto omega4_adaptator = [&beta, &ad, &weights, &grid](std::vector<double> xs) -> double {
    std::vector<double> taus = make_x_to_tau(beta)(xs);
    double weight            = 1.0;
    for (auto [i, x_i] : itertools::enumerate(xs)) {
      for (auto [j, x_j] : itertools::enumerate(grid[i])) {
        if (x_i == x_j) {
          weight *= weights[j];
          break;
        }
      }
    }
    double res = dimer_Omega4(ad, beta, taus) * weight;
    if (std::isnan(res)) {
      std::cout << "NaN encountered for taus: ";
      for (auto t : taus) { std::cout << t << " "; }
      std::cout << std::endl;
      return 0.0;
    }
    return res;
  };

  xfac::TensorCI2Param params;
  params.bondDim = bond_dim;
  xfac::CTensorCI2<double, double> TCI2_func(omega4_adaptator, grid, params);

  int counter    = 0;
  int stop_steps = 5;
  int n_iter     = 10;
  nda::array<int, 1> bond_dims(order - 1, 0);
  for (int iter = 0; iter < n_iter; iter++) {

    if (counter >= stop_steps) { break; }
    TCI2_func.iterate();
    bool flag = true;
    for (int i = 0; i < order - 1; i++) {

      auto pivots = TCI2_func.getPivotsAt(i);
      flag &= (bond_dims[i] == pivots.size());
      bond_dims[i] = pivots.size();
    }

    std::cout << "Iteration " << iter << ", bond dimension: " << bond_dims << " Error after iteration " << TCI2_func.pivotError.back() << std::endl;
    counter = flag ? counter + 1 : 0;
  }

  double sum_result = sum_TCI(TCI2_func, grid) * std::pow(beta, order);
  std::cout << "Result from TCI2: " << sum_result << std::endl;

  double exact = dimer_Omega2(ad, beta, beta / 500.0) * beta;
  std::cout << "Exact result: " << exact << std::endl;
}

// auto func = [](std::vector<double> x) -> double { return 1 + log(1 + x[0] * x[1]); };

// for (auto index1 : grid[0]) {
//   for (auto index2 : grid[1]) {
//     auto val = func({index1, index2});
//     std::cout << val << " ";
//   }
//   std::cout << std::endl;
// }

// xfac::TensorCI2Param params;
// params.bondDim = bond_dim;
// xfac::CTensorCI2<double, double> TCI2_func(func, grid, params);
// int n_iter = 10;
// int n_piv  = 0;
// for (int iter = 0; iter < n_iter; iter++) {
//   TCI2_func.iterate();
//   std::cout << "Iteration " << iter << " completed." << std::endl;
//   auto pivots = TCI2_func.getPivotsAt(0);

//   for (auto pivot : pivots) {
//     for (int j = 0; j < pivot.size(); j++) { std::cout << pivot[j] << " "; }
//     std::cout << std::endl;
//   }
//   if (n_piv == pivots.size()) { break; }
//   n_piv = pivots.size();
// }

// arma::mat integ2(n_piv, 1, arma::fill::zeros);
// // arma::mat P(2, 2);
// // P(0, 0) = 1.64;
// // P(0, 1) = 1;
// // P(1, 0) = 1;
// // P(1, 1) = 1;

// arma::mat integ(1, n_piv, arma::fill::zeros);
// // std::cout << "T0" << std::endl;
// // for (auto [index,w] : itertools::zip(grid[0], weights)) {
// //   auto T = TCI2_func.get_T_at(0, index);
// //   std::cout << T;
// // }

// std::cout << "TP1*P" << std::endl;
// for (auto [index, w] : itertools::zip(grid[0], weights)) {
//   auto TP1 = TCI2_func.get_TP1_at(0, index);
//   // std::cout << TP1 * P;
//   integ += TP1 * w;
// }

// // std::cout << "Integ 1: " << integ << std::endl;

// //std::cout << "T1.T" << std::endl;
// for (auto [index, w] : itertools::zip(grid[1], weights)) {
//   auto T1 = TCI2_func.get_T_at(1, index);
//   //std::cout << T1.t();
//   integ2 += T1 * w;
// };

// //WARNING P1T does not work!

// // std::cout << "P*P1T" << std::endl;
// // for (auto index : grid[1]) {
// //   auto P1T = TCI2_func.get_P1T_at(1, index);
// //   std::cout << (P * P1T).t();
// // };

// std::cout << "Integ 2: " << integ2 << std::endl;
// std::cout << "Integral: " << integ * integ2 << std::endl;

// // for (auto index1 : grid[0]) {
// //   for (auto index2 : grid[1]) {
// //     auto val = TCI2_func.get_CTensorTrain().eval({index1, index2});
// //     std::cout << val << " ";
// //   }
// //   std::cout << std::endl;
// // }
