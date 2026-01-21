#include "sc_expansion/configuration.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
#include "sc_expansion/dimer_order_4.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <iostream>
#include <triqs/utility/callbacks.hpp>
#include <triqs/stat/accumulator.hpp>
#include <triqs/stat/jackknife.hpp>
#include <chrono>

std::vector<sc_expansion::adjmat> diagram_mats_2() {

  sc_expansion::adjmat D2a = {{0, 1}, {1, 0}}; //2-cycle
  return {D2a};
}

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

int main(int argc, char *argv[]) {

  int n_cycles = std::stoi(argv[1]);

  // initialize mpi
  mpi::environment env(argc, argv);
  mpi::communicator world;

  // greeting
  if (world.rank() == 0) std::cout << "Strong Coupling Monte Carlo" << std::endl;

  auto start_time = std::chrono::high_resolution_clock::now();

  // Prepare the MC parameters
  //int n_cycles            = 50000;
  int length_cycle        = 1;
  int n_warmup_cycles     = 10000;
  std::string random_name = "";
  int random_seed         = 321865 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);
  int n_bins              = 50; // Initial default, will be overridden for analysis

  //diagram mats
  auto mats2 = diagram_mats_2();
  auto mats4 = diagram_mats_4();
  auto mats6 = diagram_mats_6();

  // Construct a Monte Carlo loop
  triqs::mc_tools::mc_generic<double> StrongCouplingMC(random_name, random_seed, verbosity);

  // parameters of the model
  double U    = 8.0;
  double beta = 1.0;
  double mu   = 2.0;
  int order   = 6;

  // double exact_result = -5.57770720e-03;

  double exact_result = 2.74033646e-04;

  // construct configuration
  Configuration config(U, beta, mu, mats6, order);

  //construct the reference weight function
  auto o2 = sc_expansion::order2(U, mu, beta);
  auto o4 = sc_expansion::order4(U, mu, beta);

  //reference wieight for 4th order: o2*o2
  // auto compute_reference_weight = [&o2](std::vector<double> const &taus) {
  //   return o2.compute_sum_diagrams({taus[0], taus[1]}) * o2.compute_sum_diagrams({taus[2], taus[3]});
  // };
  // double reference_integral = -8.48378682e-02 * -8.48378682e-02;

  //reference weight for 6th order: o2*o4
  auto compute_reference_weight = [&o2, &o4](std::vector<double> const &taus) {
    return -o2.compute_sum_diagrams({taus[0], taus[1]}) * -o4.compute_sum_diagrams({taus[2], taus[3], taus[4], taus[5]});
  };
  double reference_integral = (-8.48378682e-02) * (-5.57770720e-03);

  config.set_reference_weight_function(compute_reference_weight);

  // Storage for measure
  std::vector<double> signs_vec;
  std::vector<double> refs_vec;
  signs_vec.reserve(n_cycles);
  refs_vec.reserve(n_cycles);
  nda::array<double, 1> signs;
  nda::array<double, 1> reference_vals;

  //add moves and measures
  measure my_measure(&config, signs_vec, refs_vec, signs, reference_vals);
  StrongCouplingMC.add_move(move(&config, StrongCouplingMC.get_rng()), "spin flip");
  StrongCouplingMC.add_measure(my_measure, "measure sign");

  // Run and collect results
  StrongCouplingMC.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(-1));

  if (world.rank() == 0) std::cout << "Collecting results..." << std::endl;
  StrongCouplingMC.collect_results(world);

  if (world.rank() == 0) {
    std::cout << "Chain completed successfully." << std::endl;
    std::cout << "Exact Result: " << exact_result << std::endl;

    // Dynamic binning logic
    long total_samples = signs.size();
    long block_size    = 2000;                                          // Heuristic for correlation time
    if (total_samples < block_size * 2) block_size = total_samples / 5; // Fallback for small runs
    if (block_size < 1) block_size = 1;

    long calculated_n_bins = total_samples / block_size;

    std::cout << "Analysis: Total Samples = " << total_samples << ", Block Size = " << block_size << ", Bins = " << calculated_n_bins << std::endl;

    // Use shared storage to populate accumulators
    // capacity_per_bin must be sufficient to hold block_size
    triqs::stat::accumulator<double> acc_sign(0.0, 0, calculated_n_bins, block_size + 100);
    triqs::stat::accumulator<double> acc_ref(0.0, 0, calculated_n_bins, block_size + 100);

    // Access the shared nda::arrays
    for (auto s : signs) acc_sign << s;
    for (auto r : reference_vals) acc_ref << r;

    double avg_sign = 0.0;
    for (auto s : signs) avg_sign += s;
    avg_sign /= signs.size();

    double avg_ref = 0.0;
    for (auto r : reference_vals) avg_ref += r;
    avg_ref /= reference_vals.size();

    std::cout << "Explicit Ratio Result: " << reference_integral * avg_sign / avg_ref << std::endl;

    auto result = triqs::stat::jackknife(
       [reference_integral](double avg_sign, double measured_ref) { return avg_sign * reference_integral / measured_ref; }, acc_sign, acc_ref);

    std::cout << "Jackknife mean: " << std::get<0>(result) << std::endl;
    std::cout << "Jackknife error: " << std::get<1>(result) << std::endl;

    auto end_time                         = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    double total_time                     = elapsed.count();
    long total_steps                      = n_warmup_cycles + n_cycles;
    double time_per_step                  = total_time / total_steps;

    std::cout << "Total runtime: " << total_time << " s" << std::endl;
    std::cout << "Total MC steps (warmup + cycles): " << total_steps << std::endl;
    std::cout << "Time per step: " << time_per_step << " s" << std::endl;
  }
}