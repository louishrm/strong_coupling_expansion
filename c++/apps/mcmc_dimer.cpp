#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/diagmc_configuration.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure_dimer.hpp"
#include "sc_expansion/diagmc_move.hpp"
#include "sc_expansion/measure_diagmc.hpp"
#include "sc_expansion/dual.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <h5/h5.hpp>
#include <filesystem>
#include <iostream>
#include <chrono>
#include <memory>

// Broadcast a vector of Graph objects from rank 0 to all other ranks.
static void broadcast_graphs(std::vector<sc_expansion::Graph> &graphs, mpi::communicator &world) {

  int n_graphs = (int)graphs.size();
  MPI_Bcast(&n_graphs, 1, MPI_INT, 0, world.get());

  if (world.rank() != 0) { graphs.clear(); }

  for (int g = 0; g < n_graphs; g++) {
    int V, aut, sym, fm;
    int bip_only;
    std::vector<uint8_t> adj;

    if (world.rank() == 0) {
      auto const &graph = graphs[g];
      V                 = graph.get_V();
      aut               = graph.get_automorphism_count();
      sym               = (int)graph.get_symmetry_factor();
      fm                = (int)graph.get_free_multiplicity();
      bip_only          = graph.get_bipartite_only() ? 1 : 0;
      adj               = graph.get_canonical_form();
    }

    MPI_Bcast(&V, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&aut, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&sym, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&fm, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&bip_only, 1, MPI_INT, 0, world.get());

    int adj_size = V * V;
    if (world.rank() != 0) { adj.resize(adj_size); }
    MPI_Bcast(adj.data(), adj_size, MPI_UNSIGNED_CHAR, 0, world.get());

    if (world.rank() != 0) { graphs.emplace_back(adj, V, aut, sym, fm, bip_only != 0); }
  }
}

template <typename T>
void run(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, double t_hop, double alpha, int n_warmup_cycles,
         int length_cycle, std::string random_name, int random_seed, int verbosity, bool use_diagmc) {

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    params = {Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(t_hop, 0.0), true};
  } else {
    params = {U, beta, mu, t_hop, true};
  }

  // --- Phase 1: Rank 0 generates all vacuum diagrams, then broadcasts to all ranks ---
  std::vector<sc_expansion::Graph> graphs;
  if (world.rank() == 0) {
    std::cout << "Generating vacuum diagrams on rank 0..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    sc_expansion::VacuumDiagramGenerator gen(order, params.bipartite);
    gen.generate();
    graphs  = gen.get_unique_graphs();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Generated " << graphs.size() << " unique diagrams in " << std::chrono::duration<double>(t1 - t0).count() << " s." << std::endl;
  }
  broadcast_graphs(graphs, world);

  // --- Phase 2: All ranks construct the calculator from pre-built graphs ---
  sc_expansion::FreeEnergyCalculator<2, T> calculator(params, order, graphs);

  // Print diagram structure diagnostics
  if (world.rank() == 0) {
    auto const &diags = calculator.get_diagrams();
    auto const &gr    = calculator.get_graphs();
    std::cout << "\n--- Diagram structure ---" << std::endl;
    std::cout << std::left;
    std::cout << "Graph   V     SymFactor   SpatCfgs    GlobalCfgs  TotalWeight   LocalStates" << std::endl;
    std::cout << std::string(90, '-') << std::endl;
    for (size_t d = 0; d < diags.size(); d++) {
      auto const &g  = gr[d];
      auto const &di = diags[d];
      double tw      = 0;
      for (auto const &c : di.get_valid_configurations()) tw += c.weight;
      auto ls_counts = di.get_local_state_counts();
      std::cout << d << "\t" << g.get_V() << "\t" << g.get_symmetry_factor() << "\t\t" << di.get_spatial_configurations().size() << "\t\t"
                << di.get_valid_configurations().size() << "\t\t" << tw << "\t\t[";
      for (size_t v = 0; v < ls_counts.size(); v++) { std::cout << ls_counts[v] << (v + 1 < ls_counts.size() ? "," : ""); }
      std::cout << "]" << std::endl;
    }
    std::cout << std::endl;
  }

  // --- Phase 3: MC sampling ---

  // Ensure results directory exists
  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  int n_bins     = 50;
  int block_size = (n_cycles / n_bins) + 1;

  // Shared result pointer — filled by whichever MC scheme runs
  std::shared_ptr<DimerMeasureResult> result;
  double config_U, config_beta;

  if (use_diagmc) {
    // --- DiagMC: sample over (diagram, tau) jointly ---
    auto config = std::make_unique<DiagMCConfiguration<T>>(params, order, calculator, alpha);
    config_U    = config->get_U();
    config_beta = config->beta;

    if (world.rank() == 0) { std::cout << "\nUsing DiagMC mode with " << config->get_n_diagrams() << " physical diagrams" << std::endl; }

    triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);
    MeasureDiagMC<T> meas(config.get(), n_bins, block_size, verbosity);
    mc.add_move(DiagMCMove<T>(config.get(), mc.get_rng()), "diagmc_move");
    mc.add_measure(meas, "diagmc_measure");

    auto start_time = std::chrono::high_resolution_clock::now();
    mc.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
    mc.collect_results(world);
    result = meas.result;

    if (world.rank() == 0) {
      auto end_time                         = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = end_time - start_time;
      double total_time                     = elapsed.count();
      long total_steps                      = (long)(n_warmup_cycles + n_cycles) * world.size();
      std::cout << "Total time (s): " << total_time << std::endl;
      std::cout << "Time per step (s): " << total_time / total_steps << std::endl;
    }

  } else {
    // --- Standard defensive scheme: evaluate all diagrams at every step ---
    auto config = std::make_unique<DimerConfiguration<T>>(params, order, calculator, alpha);
    config_U    = config->get_U();
    config_beta = config->beta;

    triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);
    measure_dimer<T> meas(config.get(), n_bins, block_size, mu, verbosity);
    mc.add_move(move<T>(config.get(), mc.get_rng()), "time_swap");
    mc.add_measure(meas, "dimer_measure");

    auto start_time = std::chrono::high_resolution_clock::now();
    mc.warmup_and_accumulate(n_warmup_cycles, n_cycles, length_cycle, triqs::utility::clock_callback(-1));
    mc.collect_results(world);
    result = meas.result;

    if (world.rank() == 0) {
      auto end_time                         = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> elapsed = end_time - start_time;
      double total_time                     = elapsed.count();
      long total_steps                      = (long)(n_warmup_cycles + n_cycles) * world.size();
      std::cout << "Total time (s): " << total_time << std::endl;
      std::cout << "Time per step (s): " << total_time / total_steps << std::endl;
    }
  }

  if (world.rank() == 0) {
    // Save results to HDF5
    std::string filename = "./results/dimer_data_order_" + std::to_string(order) + "_U_" + std::to_string(config_U) + "_beta_"
       + std::to_string(config_beta) + "_mu_" + std::to_string(mu) + ".h5";
    h5::file file(filename, 'w');
    h5_write(file, "mean", result->coeff);
    h5_write(file, "error", result->error);
    h5_write(file, "mean_sign", result->mean_sign);
    h5_write(file, "sign_error", result->sign_error);
    h5_write(file, "mean_abs_integrand", result->mean_abs);
    h5_write(file, "abs_integrand_error", result->abs_error);
    h5_write(file, "mu", mu);
  }
}

int main(int argc, char *argv[]) {

  if (argc < 8) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu t_hop alpha [use_dual] [use_diagmc]" << std::endl;
    }
    return 1;
  }

  int order       = std::stoi(argv[1]);
  int n_cycles    = std::stoi(argv[2]);
  double U        = std::stod(argv[3]);
  double beta     = std::stod(argv[4]);
  double mu       = std::stod(argv[5]);
  double t_hop    = std::stod(argv[6]);
  double alpha    = std::stod(argv[7]);
  bool use_dual   = (argc > 8 ? std::stoi(argv[8]) != 0 : false);
  bool use_diagmc = (argc > 9 ? std::stoi(argv[9]) != 0 : false);

  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    std::cout << "=== Strong Coupling MC (Dimer) ===" << std::endl;
    std::cout << "MPI ranks: " << world.size() << std::endl;
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " t_hop=" << t_hop << " alpha=" << alpha
              << " diagmc=" << use_diagmc << std::endl;
  }

  int length_cycle        = 1;
  int n_warmup_cycles     = 5000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (use_dual) {
    run<Dual>(world, order, n_cycles, U, beta, mu, t_hop, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity, use_diagmc);
  } else {
    run<double>(world, order, n_cycles, U, beta, mu, t_hop, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity, use_diagmc);
  }

  return 0;
}
