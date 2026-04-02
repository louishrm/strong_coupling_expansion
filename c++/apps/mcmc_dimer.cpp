#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include "sc_expansion/dual.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <iostream>
#include <chrono>
#include <memory>

template <typename T>
void run(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, double t_hop, int n_warmup_cycles, int length_cycle,
         std::string random_name, int random_seed, int verbosity) {

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    params = {Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(t_hop, 0.0), true};
  } else {
    params = {U, beta, mu, t_hop, true};
  }

  // Construct dimer configuration — no reference integral needed
  auto config = std::make_unique<DimerConfiguration<T>>(params, order);

  // Construct MC loop
  triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  int measure_seed = random_seed + 99871234;
  measure_dimer<T> meas(config.get(), n_bins, block_size, mu, measure_seed);
  mc.add_move(move<T>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_measure");

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
  }
}

int main(int argc, char *argv[]) {

  if (argc < 7) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu t_hop [use_dual]" << std::endl;
    }
    return 1;
  }

  int order     = std::stoi(argv[1]);
  int n_cycles  = std::stoi(argv[2]);
  double U      = std::stod(argv[3]);
  double beta   = std::stod(argv[4]);
  double mu     = std::stod(argv[5]);
  double t_hop  = std::stod(argv[6]);
  bool use_dual = (argc > 7 ? std::stoi(argv[7]) != 0 : false);

  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    std::cout << "=== Strong Coupling MC (Dimer) ===" << std::endl;
    std::cout << "MPI ranks: " << world.size() << std::endl;
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " t_hop=" << t_hop << std::endl;
  }

  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (use_dual) {
    run<Dual>(world, order, n_cycles, U, beta, mu, t_hop, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  } else {
    run<double>(world, order, n_cycles, U, beta, mu, t_hop, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  }

  return 0;
}
