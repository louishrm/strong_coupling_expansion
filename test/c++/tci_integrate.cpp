#include <cmath>
#include <iostream>
#include <xfac/tensor/tensor_ci_2.h>
#include <nda/nda.hpp>
#include "../../c++/sc_expansion/diagram.hpp"
#include "../../c++/sc_expansion/hubbard_atom.hpp"
#include "../../c++/sc_expansion/cumulant.hpp"

double dimer_Omega_n(auto ad, double beta, std::vector<double> tau, std::vector<Diagram> Ds) {

  double res = 0.0;
  for (const auto &D : Ds) { res += D.evaluate_at_taus(ad, beta, tau); }
  return res;
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
  double beta = 6.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops     = hubbard_atom::make_fops();
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
  for (int i = 0; i < grid_size; i++) { weights.push_back(1.0 / (grid_size)); }

  //Diagram D({{0, 1}, {1, 0}});

  adjmat D4a                    = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  adjmat D4b                    = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
  adjmat D4c                    = {{0, 2}, {2, 0}};                                         //2-cycle with double lines
  std::vector<Diagram> diagrams = {Diagram(D4a), Diagram(D4b), Diagram(D4c)};

  //std::vector<Diagram> diagrams = {Diagram({{0, 1}, {1, 0}})}; // second order diagram

  auto dimer_omega_n_adaptator = [&beta, &ad, &weights, &grid, &diagrams](std::vector<double> xs) -> double {
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
    double res = dimer_Omega_n(ad, beta, taus, diagrams) * weight;
    return res;
  };

  xfac::TensorCI2Param params;
  params.bondDim = bond_dim;
  xfac::CTensorCI2<double, double> TCI2_func(dimer_omega_n_adaptator, grid, params);

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
  return 0;
}
