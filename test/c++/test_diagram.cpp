#include "../../c++/sc_expansion/diagram.hpp"
#include "../../c++/sc_expansion/hubbard_atom.hpp"
#include "../../c++/sc_expansion/cumulant.hpp"

#include <fstream>
#include <iostream>
#include <nda/nda.hpp>

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

double dimer_Omega2_new(auto ad, double beta, Diagram &D, double delta) {

  double result = 0.0;

  for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {

    std::vector<double> tau = {0, tau2};
    double result_a         = D.evaluate_at_taus(ad, beta, tau);
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

  double C02_up = compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 0}}, {{tau[2], 0}, {tau[3], 0}}, ad, beta); //C(2 4 | 1 3)

  double G01_14_up = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[0], 0}}, ad, beta);
  double G01_32_up = compute_cumulant_decomposition({{tau[3], 0}}, {{tau[1], 0}}, ad, beta);
  result += 2 * C02_up * G01_14_up * G01_32_up;

  double C02_updown = compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 1}}, {{tau[2], 0}, {tau[3], 1}}, ad, beta); //C(2 4 | 1 3)

  double G01_14_up_2 = compute_cumulant_decomposition({{tau[2], 0}}, {{tau[0], 0}}, ad, beta); //G(1|4)
  double G01_32_down = compute_cumulant_decomposition({{tau[1], 1}}, {{tau[3], 1}}, ad, beta); //G(3|2)

  result += 2 * C02_updown * G01_14_up_2 * G01_32_down;

  double sign              = 1.0;
  double symmetry_factor   = 1.0 / 2.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega4c(auto ad, double beta, std::vector<double> tau) {

  double result = 0.0;

  result += 2 * compute_cumulant_decomposition({{tau[0], 0}, {tau[2], 0}}, {{tau[2], 0}, {tau[3], 0}}, ad, beta)
     * compute_cumulant_decomposition({{tau[2], 0}, {tau[3], 0}}, {{tau[0], 0}, {tau[1], 0}}, ad, beta);

  result += 2 * compute_cumulant_decomposition({{tau[0], 0}, {tau[1], 1}}, {{tau[2], 0}, {tau[3], 1}}, ad, beta)
     * compute_cumulant_decomposition({{tau[2], 0}, {tau[3], 1}}, {{tau[0], 0}, {tau[1], 1}}, ad, beta);

  result += 2 * compute_cumulant_decomposition({{tau[0], 1}, {tau[1], 0}}, {{tau[2], 0}, {tau[3], 1}}, ad, beta)
     * compute_cumulant_decomposition({{tau[2], 0}, {tau[3], 1}}, {{tau[0], 1}, {tau[1], 0}}, ad, beta);

  double sign              = 1.0;
  double symmetry_factor   = 1 / 8.0;
  double free_multiplicity = 2.0;
  double prefactor         = sign * symmetry_factor * free_multiplicity;

  return prefactor * result;
}

double dimer_Omega4(auto ad, double beta, double delta) {

  double result = 0.0;

  for (double tau2 = 0.0; tau2 <= beta; tau2 += delta) {
    for (double tau3 = 0.0; tau3 <= beta; tau3 += delta) {
      for (double tau4 = 0.0; tau4 <= beta; tau4 += delta) {

        std::vector<double> tau = {0, tau2, tau3, tau4};

        double result_a = dimer_Omega4a(ad, beta, tau);
        double result_b = dimer_Omega4b(ad, beta, tau);
        double result_c = dimer_Omega4c(ad, beta, tau);
        result += (result_a + result_b + result_c) * delta * delta * delta;
      }
    }
  }
  return result;
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

template <typename T> inline void save(const std::string filename, const std::string name, const T &object) {
  h5::file hfile(filename, 'a');
  h5::write(hfile, name, object);
  hfile.close();
}

// 1. Define the comparator
struct NearLess {
  double epsilon = 1e-6; // Adjust this threshold based on your needs
  bool operator()(const double &a, const double &b) const {
    // Only return true if a is definitely smaller than b (outside epsilon)
    return (b - a) > epsilon;
  }
};

int main() {

  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;

  triqs::hilbert_space::fundamental_operator_set fops     = hubbard_atom::make_fops();
  triqs::operators::many_body_operator_generic<double> H0 = hubbard_atom::make_H0(U, mu);
  triqs::atom_diag::atom_diag<false> ad(H0, fops, {}); // atom_diag object

  Diagram D({{0, 1}, {1, 0}}, U, beta, mu); // second order diagram

  int grid_size = 20;
  int order     = 4;
  std::vector<std::vector<double>> grid;

  for (int i = 0; i < order; i++) {
    auto row = linspace(0, beta, grid_size + 2, true);
    row.pop_back();         // remove the endpoint 1.0
    row.erase(row.begin()); // remove the start point 0.0
    grid.push_back(row);
  }

  // std::string directory = "/Users/louissharma/Desktop/tci_python/";
  // std::string name      = "dimer_Omega2a_U" + std::to_string(int(U)) + "_mu" + std::to_string(int(mu)) + "_beta" + std::to_string(int(beta))
  //    + std::to_string("_tdiff") + "_size" + std::to_string(grid_size);

  // std::ofstream outfile(directory + name);
  // //double exact_integ = 0.0;
  // nda::matrix<double> omega2a_mat(grid_size, grid_size);
  // for (size_t i = 0; i < grid[0].size(); ++i) {
  //   double x1 = grid[0][i];
  //   for (size_t j = 0; j < grid[1].size(); ++j) {
  //     double x2 = grid[1][j];

  //     double tau1 = x1 + x2;
  //     double tau2 = x2;

  //     // std::vector<double> x = {x1, x2};
  //     // auto taus             = make_x_to_tau(beta)(x);
  //     double val = D.evaluate_at_taus(ad, beta, {tau1, tau2});

  //     omega2a_mat(i, j) = val;
  //     //exact_integ += val * 1 / 50.0 * 1 / 50.0;
  //   }
  // }
  // //std::cout << "Exact integration value for Omega2a: " << exact_integ << std::endl;
  // save(directory + name, "omega2a", omega2a_mat);

  // outfile.close();

  // double F2_new = dimer_Omega2_new(ad, beta, D, beta / 1000.0);
  // std::cout << "Free energy 2nd order from Diagram class for U: " << U << ", mu: " << mu << ", beta: " << beta << " is " << F2_new * beta
  //           << std::endl;

  // double delta = beta / 1000.0;
  // double F2    = dimer_Omega2(ad, beta, delta);
  // std::cout << "Free energy 2nd order for U: " << U << ", mu: " << mu << ", beta: " << beta << " is " << F2 * beta << std::endl;

  adjmat D4a                    = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  adjmat D4b                    = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
  adjmat D4c                    = {{0, 2}, {2, 0}};                                         //2-cycle with double lines
  std::vector<Diagram> diagrams = {Diagram(D4a, U, beta, mu), Diagram(D4b, U, beta, mu), Diagram(D4c, U, beta, mu)};

  // std::string name = "dimer_Omega4_U" + std::to_string(int(U)) + "_mu" + std::to_string(int(mu)) + "_beta" + std::to_string(int(beta)) + "_size"
  //    + std::to_string(grid_size) + ".txt";

  // std::ofstream outfile(directory + name);

  // outfile.precision(16);

  // nda::array<double, 4> omega4_arr(grid_size, grid_size, grid_size, grid_size);

  // for (size_t i = 0; i < grid[0].size(); ++i) {
  //   double x1 = grid[0][i];

  //   // j index for tau2
  //   for (size_t j = 0; j < grid[1].size(); ++j) {
  //     double x2 = grid[1][j];

  //     // k index for tau3
  //     for (size_t k = 0; k < grid[2].size(); ++k) {
  //       double x3 = grid[2][k];

  //       for (size_t l = 0; l < grid[3].size(); l++) {

  //         double x4 = grid[3][l];

  //         double tau1              = x1 + x2 + x3 + x4;
  //         double tau2              = x2 + x3 + x4;
  //         double tau3              = x3 + x4;
  //         double tau4              = x4;
  //         std::vector<double> taus = {tau1, tau2, tau3, tau4};
  //         double sum_val           = 0.0;
  //         for (size_t diag_idx = 0; diag_idx < diagrams.size(); diag_idx++) {
  //           // Evaluate each diagram contribution
  //           double val = diagrams[diag_idx].evaluate_at_taus(ad, beta, taus);
  //           sum_val += val; // Sum the contributions
  //         }
  //         omega4_arr(i, j, k, l) = sum_val;
  //       } //end l
  //     } // end k
  //   } // end j
  // } // end i

  // save(directory + name, "omega4", omega4_arr);
  // outfile.close();

  // std::string name = "dimer_Omega4_U_3tensor" + std::to_string(int(U)) + "_mu" + std::to_string(int(mu)) + "_beta" + std::to_string(int(beta))
  //    + "_size" + std::to_string(grid_size) + ".txt";

  // std::ofstream outfile(directory + name);

  // outfile.precision(16);

  // nda::array<double, 3> omega4_arr(grid_size, grid_size, grid_size);

  // for (size_t i = 0; i < grid[0].size(); ++i) {
  //   double x1 = grid[0][i];

  //   // j index for tau2
  //   for (size_t j = 0; j < grid[1].size(); ++j) {
  //     double x2 = grid[1][j];

  //     // k index for tau3
  //     for (size_t k = 0; k < grid[2].size(); ++k) {
  //       double x3 = grid[2][k];

  //       double tau1              = x1 + x2 + x3;
  //       double tau2              = x2 + x3;
  //       double tau3              = x3;
  //       double tau4              = 0.0;
  //       std::vector<double> taus = {tau1, tau2, tau3, tau4};
  //       double sum_val           = 0.0;
  //       for (size_t diag_idx = 0; diag_idx < diagrams.size(); diag_idx++) {
  //         // Evaluate each diagram contribution
  //         double val = diagrams[diag_idx].evaluate_at_taus(ad, beta, taus);
  //         sum_val += val; // Sum the contributions
  //       }
  //       omega4_arr(i, j, k) = sum_val;
  //     } // end k
  //   } // end j
  // } // end i

  // save(directory + name, "omega4", omega4_arr);
  // outfile.close();

  double integ = 0.0;
  double delta = beta / grid_size;
  for (size_t i = 0; i < grid[0].size(); ++i) {
    double tau1 = grid[0][i];

    // j index for tau2
    for (size_t j = 0; j < grid[1].size(); ++j) {
      double tau2 = grid[1][j];

      // k index for tau3
      for (size_t k = 0; k < grid[2].size(); ++k) {
        double tau3 = grid[2][k];

        for (size_t l = 0; l < grid[3].size(); l++) {

          double tau4              = grid[3][l];
          std::vector<double> taus = {tau1, tau2, tau3, tau4};
          double sum_val           = 0.0;
          for (size_t diag_idx = 0; diag_idx < diagrams.size(); diag_idx++) {
            // Evaluate each diagram contribution
            double val = diagrams[diag_idx].evaluate_at_taus(taus);
            sum_val += val; // Sum the contributions
          }
          integ += sum_val * delta * delta * delta * delta;
        } //end l
      } // end k
    } // end j
  } // end i

  std::cout << "Free energy 4th order from Diagram class for U: " << U << ", mu: " << mu << ", beta: " << beta << " is " << integ * beta << std::endl;

  std::vector<double> tau_set = {0.356, 0.679, 0.709, 0.9};

  std::sort(tau_set.begin(), tau_set.end());
  std::set<double, NearLess> unique_vals;
  do {
    double val = diagrams[2].evaluate_at_taus(tau_set);
    for (size_t i = 0; i < tau_set.size(); i++) { std::cout << tau_set[i] << ", "; }
    std::cout << "}: " << val << std::endl;
    unique_vals.insert(val);
  } while (std::next_permutation(tau_set.begin(), tau_set.end()));

  std::cout << "Unique values count: " << unique_vals.size() << std::endl;
  for (const auto &val : unique_vals) { std::cout << "Unique value: " << val << std::endl; }
  return 0;
}