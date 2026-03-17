#include "sc_expansion/configuration.hpp"
#include "sc_expansion/free_energy_order.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include "sc_expansion/dual.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <iostream>
#include <triqs/utility/callbacks.hpp>
#include <triqs/stat/accumulator.hpp>
#include <triqs/stat/jackknife.hpp>
#include <h5/h5.hpp>
#include <chrono>
#include <numeric>
#include <algorithm>

template <typename T, typename Order> std::pair<double, double> compute_exact_integral_infinite_U(Order &o, int n, double beta) {
  //Return \int |U_inf| d tau and \int U_inf d tau (both abs and signed version)
  std::vector<double> taus(n);
  std::iota(taus.begin(), taus.end(), 0.0);

  double sum_abs    = 0.0;
  double sum_signed = 0.0;
  do {
    if constexpr (std::is_same_v<T, Dual>) {
      double val = o.compute_sum_diagrams(taus, true, false).derivative;
      sum_abs += std::abs(val);
      sum_signed += val;
    } else {
      double val = o.compute_sum_diagrams(taus, true, false);
      sum_abs += std::abs(val);
      sum_signed += val;
    }
  } while (std::next_permutation(taus.begin(), taus.end()));

  double fact = 1.0;
  for (int i = 1; i <= n; ++i) fact *= i;

  std::pair<double, double> result;
  result.first  = (std::pow(beta, n) / fact) * sum_abs;
  result.second = (std::pow(beta, n) / fact) * sum_signed;
  return result;
}

template <typename T>
void run_mc(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, bool bipartite, double alpha, int n_warmup_cycles,
            int length_cycle, std::string random_name, int random_seed, int verbosity) {

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    params = {Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), bipartite};
  } else {
    params = {U, beta, mu, bipartite};
  }

  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    std::cout << "Computing reference integral on master rank..." << std::endl;
    sc_expansion::FreeEnergyCalculator<T> calculator(params, order);
    std::pair<double, double> reference_vals = compute_exact_integral_infinite_U<T>(calculator, order, beta);
    reference_integral                       = reference_vals.first;
    signed_reference_integral                = reference_vals.second;
    std::cout << "Done computing reference integral. Value: " << signed_reference_integral << std::endl;
  }

  // Broadcast reference values to all ranks
  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<double> StrongCouplingMC(random_name, random_seed, verbosity);

  // construct configuration
  Configuration<T> config(params, order, alpha);

  int n_bins     = 50; // Standard choice for Jackknife
  int block_size = (n_cycles / n_bins) + 1;

  // add moves and measures
  measure<T> my_measure(&config, reference_integral, signed_reference_integral, n_bins, block_size, mu);
  StrongCouplingMC.add_move(move<T>(&config, StrongCouplingMC.get_rng()), "time_swap");
  StrongCouplingMC.add_measure(my_measure, "defensive_measure");

  auto start_time = std::chrono::high_resolution_clock::now();
  StrongCouplingMC.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(-1));

  StrongCouplingMC.collect_results(world);

  if (world.rank() == 0) {
    auto end_time                         = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    double total_time                     = elapsed.count();
    long total_steps                      = (long)(n_warmup_cycles + n_cycles) * world.size();
    double time_per_step                  = total_time / total_steps;

    std::cout << "Total time (s): " << total_time << std::endl;
    std::cout << "Time per step (s): " << time_per_step << std::endl;
    std::cout << "Steps per second: " << 1.0 / time_per_step << std::endl;

    std::cout << "Exact result (Infinite U, Order " << order << "): " << signed_reference_integral << std::endl;
  }
}

int main(int argc, char *argv[]) {

  if (argc < 6) {
    if (mpi::communicator().rank() == 0) { std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu [bipartite] [alpha] [use_dual]" << std::endl; }
    return 1;
  }

  int order      = std::stoi(argv[1]);
  int n_cycles   = std::stoi(argv[2]);
  double U       = std::stod(argv[3]);
  double beta    = std::stod(argv[4]);
  double mu      = std::stod(argv[5]);
  bool bipartite = (argc > 6 ? std::stoi(argv[6]) != 0 : true);
  double alpha   = (argc > 7 ? std::stod(argv[7]) : 0.5);
  bool use_dual  = (argc > 8 ? std::stoi(argv[8]) != 0 : false);

  // initialize mpi
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // greeting
  if (world.rank() == 0) {
    std::cout << "Strong Coupling Monte Carlo" << std::endl;
    std::cout << "Number of MPI processes: " << world.size() << std::endl;
    std::cout << "U=" << U << " beta=" << beta << " mu=" << mu << " bipartite=" << bipartite << " alpha=" << alpha << " use_dual=" << use_dual << std::endl;
  }

  // MC parameters
  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (use_dual) {
    run_mc<Dual>(world, order, n_cycles, U, beta, mu, bipartite, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  } else {
    run_mc<double>(world, order, n_cycles, U, beta, mu, bipartite, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  }

  return 0;
}
