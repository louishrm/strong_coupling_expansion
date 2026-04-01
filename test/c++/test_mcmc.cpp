/*
 * MCMC benchmark test for the atomic strong coupling expansion.
 *
 * Runs order-4 MCMC with free multiplicities set to 1 (matching the 2-site
 * ED benchmark from analytical/benchmark_atomic_expansion.py) and checks
 * that the MC estimate converges to the exact permutation-based coefficient.
 *
 * Usage:  mpirun -np 4 ./test_mcmc
 *         (also works with 1 rank: ./test_mcmc)
 */

#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/combinatorics.hpp"
#include <mpi/mpi.hpp>
#include <triqs/mc_tools/random_generator.hpp>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cmath>

int main(int argc, char **argv) {

  mpi::environment env(argc, argv);
  mpi::communicator world;

  // =====================================================================
  // Parameters (matching the Python ED benchmark)
  // =====================================================================
  double U    = 8.0;
  double beta = 2.0;
  double mu   = 3.0;
  int order   = 4;

  double alpha   = 0.01;
  int n_samples  = 100000;
  int n_warmup   = 5000;

  sc_expansion::Parameters<double> params{U, beta, mu, 0.0, true};

  // =====================================================================
  // Create calculator with fm=1 (2-site benchmark: 1 embedding per diagram)
  // =====================================================================
  sc_expansion::FreeEnergyCalculator<1, double> calc(params, order, /*override_fm=*/1);

  // =====================================================================
  // Exact reference: infinite-U coefficient by exhaustive permutation
  // =====================================================================
  auto [ref_abs, ref_signed] = calc.compute_infinite_U_coefficient(false);

  // =====================================================================
  // "Exact" finite-U coefficient by same permutation method
  // (This uses taus = permutations of {0,1,...,n-1} and beta^n/n! scaling)
  // =====================================================================
  double perm_sum = 0.0;
  {
    std::vector<double> perm_taus(order);
    std::iota(perm_taus.begin(), perm_taus.end(), 0.0);
    do { perm_sum += calc.compute_sum_diagrams(perm_taus, false, false); } while (std::next_permutation(perm_taus.begin(), perm_taus.end()));
  }
  double exact_coeff = std::pow(beta, order) / sc_expansion::factorial(order) * perm_sum;

  // =====================================================================
  // MCMC loop with defensive importance sampling
  // W = |f - f_inf| + alpha * |f_inf|
  // Estimator: I = I_ref * <(f-f_inf)/W> / <|f_inf|/W> + I_ref_signed
  // =====================================================================
  triqs::mc_tools::random_generator rng("mt19937", 32186222 + world.rank() * 786512);

  // Initialize random taus in [0, beta]
  std::vector<double> taus(order);
  for (auto &t : taus) t = rng(beta);

  // Evaluate initial state
  double f_val = calc.compute_sum_diagrams(taus, false, true);
  double f_inf = calc.compute_sum_diagrams(taus, true, true);
  double W     = std::abs(f_val - f_inf) + alpha * std::abs(f_inf);

  // Accumulators
  double sum_int = 0.0, sum_ref = 0.0;
  int n_acc = 0, n_accepted = 0;

  for (int step = 0; step < n_warmup + n_samples; step++) {

    // Propose: change one random tau
    int idx       = rng(order);
    double old_tau = taus[idx];
    taus[idx]     = rng(beta);

    double new_f   = calc.compute_sum_diagrams(taus, false, true);
    double new_inf = calc.compute_sum_diagrams(taus, true, true);
    double new_W   = std::abs(new_f - new_inf) + alpha * std::abs(new_inf);

    // Metropolis
    double ratio = (W > 0.0) ? new_W / W : (new_W > 0.0 ? 1.0e100 : 1.0);
    if (rng(1.0) < std::min(1.0, ratio)) {
      f_val = new_f;
      f_inf = new_inf;
      W     = new_W;
      n_accepted++;
    } else {
      taus[idx] = old_tau;
    }

    // Accumulate after warmup
    if (step >= n_warmup && W > 0.0) {
      sum_int += (f_val - f_inf) / W;
      sum_ref += std::abs(f_inf) / W;
      n_acc++;
    }
  }

  // =====================================================================
  // MPI reduction: average local estimates across ranks
  // =====================================================================
  double local_avg_int = (n_acc > 0) ? sum_int / n_acc : 0.0;
  double local_avg_ref = (n_acc > 0) ? sum_ref / n_acc : 0.0;

  double global_avg_int = mpi::all_reduce(local_avg_int, world) / world.size();
  double global_avg_ref = mpi::all_reduce(local_avg_ref, world) / world.size();

  // Ratio estimator
  double mc_estimate = (std::abs(global_avg_ref) > 1e-18) ? (global_avg_int / global_avg_ref) * ref_abs + ref_signed : ref_signed;

  // =====================================================================
  // Report and check
  // =====================================================================
  if (world.rank() == 0) {
    double rel_err = (std::abs(exact_coeff) > 1e-18) ? std::abs(mc_estimate - exact_coeff) / std::abs(exact_coeff) : std::abs(mc_estimate);

    std::cout << "=== MCMC Atomic Benchmark (order " << order << ", fm=1) ===" << std::endl;
    std::cout << "MPI ranks:          " << world.size() << std::endl;
    std::cout << "Samples per rank:   " << n_samples << std::endl;
    std::cout << "Acceptance rate:    " << (double)n_accepted / (n_warmup + n_samples) << std::endl;
    std::cout << "Exact (permutation):" << exact_coeff << std::endl;
    std::cout << "Infinite-U ref:     " << ref_signed << std::endl;
    std::cout << "MC estimate:        " << mc_estimate << std::endl;
    std::cout << "Relative error:     " << rel_err << std::endl;

    if (rel_err < 0.15) {
      std::cout << "PASS" << std::endl;
    } else {
      std::cerr << "FAIL: relative error " << rel_err << " exceeds 15% tolerance" << std::endl;
      return 1;
    }
  }

  return 0;
}
