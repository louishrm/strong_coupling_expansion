#include "sc_expansion/configuration.hpp"
#include "sc_expansion/free_energy_order.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <iostream>
#include <triqs/utility/callbacks.hpp>
#include <triqs/stat/accumulator.hpp>
#include <triqs/stat/jackknife.hpp>
#include <h5/h5.hpp>
#include <chrono>
#include <numeric>
#include <algorithm>

double Zat(double U, double beta, double mu) { return 1 + 2.0 * std::exp(beta * mu) + std::exp(beta * (2.0 * mu - U)); }

template <typename Order> double compute_exact_integral_infinite_U(Order &o, int n, double beta) {
  std::vector<double> taus(n);
  std::iota(taus.begin(), taus.end(), 0.0);

  double sum = 0.0;
  do { sum += o.compute_sum_diagrams(taus, true); } while (std::next_permutation(taus.begin(), taus.end()));

  double fact = 1.0;
  for (int i = 1; i <= n; ++i) fact *= i;

  return (std::pow(beta, n) / fact) * sum;
}

int main(int argc, char *argv[]) {

  if (argc < 6) {
    if (mpi::communicator().rank() == 0) { std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu [alpha]" << std::endl; }
    return 1;
  }

  int order    = std::stoi(argv[1]);
  int n_cycles = std::stoi(argv[2]);
  double U     = std::stod(argv[3]);
  double beta  = std::stod(argv[4]);
  double mu    = std::stod(argv[5]);
  double alpha = (argc > 6 ? std::stod(argv[6]) : 0.5);

  // initialize mpi
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // greeting
  if (world.rank() == 0) {
    std::cout << "Strong Coupling Monte Carlo" << std::endl;
    std::cout << "Number of MPI processes: " << world.size() << std::endl;
    std::cout << "U=" << U << " beta=" << beta << " mu=" << mu << " alpha=" << alpha << std::endl;
  }

  auto start_time = std::chrono::high_resolution_clock::now();

  // MC parameters
  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<double> StrongCouplingMC(random_name, random_seed, verbosity);

  // construct configuration
  sc_expansion::Parameters params{U, beta, mu};
  Configuration config(params, order, alpha);

  // parameters of the model
  double reference_integral = 0;
  sc_expansion::FreeEnergyCalculator calculator(params, order);
  reference_integral = compute_exact_integral_infinite_U(calculator, order, beta);

  long total_cycles = n_cycles * world.size(); // Total samples across all cores
  int n_bins        = 50;                      // Standard choice for Jackknife
  int block_size    = (n_cycles / n_bins) + 1;

  // add moves and measures
  measure my_measure(&config, reference_integral, n_bins, block_size, mu);
  StrongCouplingMC.add_move(move(&config, StrongCouplingMC.get_rng()), "time_swap");
  StrongCouplingMC.add_measure(my_measure, "defensive_measure");

  // Run and collect results
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

    std::cout << "Exact result (Infinite U, Order " << order << "): " << reference_integral << std::endl;
  }

  return 0;
}

/*
std::vector<sc_expansion::adjmat> diagram_mats_4() {
  sc_expansion::adjmat D4a = {{0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}}; //4-cycle
  sc_expansion::adjmat D4b = {{0, 1, 1}, {1, 0, 0}, {1, 0, 0}};                        //3-cycle with double lines
  sc_expansion::adjmat D4c = {{0, 2}, {2, 0}};                                         //2-cycle with double lines
  return {D4a, D4b, D4c};
}

std::vector<sc_expansion::adjmat> diagram_mats_6() {

  sc_expansion::adjmat D6a = {{0, 1, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0}, {0, 0, 0, 1, 0, 0},
                              {0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 1}, {1, 0, 0, 0, 0, 0}};                          //6-cycle
  sc_expansion::adjmat D6b = {{0, 3}, {3, 0}};                                                                      //watermelon triple
  sc_expansion::adjmat D6c = {{0, 1, 1, 1}, {1, 0, 0, 0}, {1, 0, 0, 0}, {1, 0, 0, 0}};                              //petal with 4 vertice
  sc_expansion::adjmat D6d = {{0, 1, 1, 0, 0}, {1, 0, 0, 0, 0}, {0, 0, 0, 1, 0}, {0, 0, 0, 0, 1}, {1, 0, 0, 0, 0}}; //square +digon
  sc_expansion::adjmat D6e = {{0, 1, 1, 0}, {1, 0, 0, 1}, {1, 0, 0, 0}, {0, 1, 0, 0}};                              //crab diagram
  sc_expansion::adjmat D6f = {{0, 2, 1}, {2, 0, 0}, {1, 0, 0}};                                                     //watermelon double + digon
  sc_expansion::adjmat D6g = {{0, 2, 0, 0}, {1, 0, 1, 0}, {0, 0, 0, 1}, {1, 0, 0, 0}};                              //square with one double line
  return {D6a, D6b, D6c, D6d, D6e, D6f, D6g};
}
*/
