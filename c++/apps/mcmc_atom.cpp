#include "sc_expansion/atomic_configuration.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include "sc_expansion/dual.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <iostream>
#include <chrono>
#include <memory>

template <typename T>
void run(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, bool bipartite, double alpha, int n_warmup_cycles,
         int length_cycle, std::string random_name, int random_seed, int verbosity) {

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    params = {Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(0.0, 0.0), bipartite};
  } else {
    params = {U, beta, mu, 0.0, bipartite};
  }

  // Compute exact infinite-U reference integral on master rank
  double reference_integral        = 0.0;
  double signed_reference_integral = 0.0;

  if (world.rank() == 0) {
    std::cout << "Computing reference integral on master rank..." << std::endl;
    sc_expansion::FreeEnergyCalculator<1, T> calculator(params, order);
    auto [ref_abs, ref_signed] = calculator.compute_infinite_U_coefficient(false);
    reference_integral         = ref_abs;
    signed_reference_integral  = ref_signed;
    std::cout << "Reference integral: " << signed_reference_integral << " (abs: " << reference_integral << ")" << std::endl;
  }

  mpi::broadcast(reference_integral, world);
  mpi::broadcast(signed_reference_integral, world);

  // Construct configuration
  auto config = std::make_unique<AtomicConfiguration<T>>(params, order, alpha);

  // Construct MC loop
  triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  measure<T> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu);
  mc.add_move(move<T>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "defensive_measure");

  auto start_time = std::chrono::high_resolution_clock::now();
  mc.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
  mc.collect_results(world);

  if (world.rank() == 0) {
    auto end_time                         = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    double total_time                     = elapsed.count();
    long total_steps                      = (long)(n_warmup_cycles + n_cycles) * world.size();

    std::cout << "Total time (s): " << total_time << std::endl;
    std::cout << "Time per step (s): " << total_time / total_steps << std::endl;
    std::cout << "Exact infinite-U coefficient (order " << order << "): " << signed_reference_integral << std::endl;
  }
}

int main(int argc, char *argv[]) {

  if (argc < 6) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu [bipartite] [alpha] [use_dual]" << std::endl;
    }
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

  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    std::cout << "=== Strong Coupling MC (Atomic) ===" << std::endl;
    std::cout << "MPI ranks: " << world.size() << std::endl;
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " bipartite=" << bipartite << " alpha=" << alpha << std::endl;
  }

  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (use_dual) {
    run<Dual>(world, order, n_cycles, U, beta, mu, bipartite, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  } else {
    run<double>(world, order, n_cycles, U, beta, mu, bipartite, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  }

  return 0;
}
