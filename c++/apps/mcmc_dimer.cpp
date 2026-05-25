// MCMC for the staggered (triangular) dimer expansion.
//
// Computes the order-n coefficient of the strong-coupling free-energy series
// for the two-site (dimer) cluster on the staggered superlattice. Every NN
// bond crosses the A/B sublattices, so the expansion runs over the FULL set
// of vacuum diagrams (bipartite + non-bipartite — odd cycles included). This
// app forces `params.bipartite = false` accordingly.
//
// Reference scheme: Configuration samples with weight W = |Omega + alpha|
// and the coefficient is recovered via the ratio estimator
//   coeff = alpha * beta^n * <Omega/W> / <alpha/W>.

#include "sc_expansion/dimer/configuration.hpp"
#include "sc_expansion/dimer/free_energy_calculator.hpp"
#include "sc_expansion/dimer/measure_dimer.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/dual.hpp"
#include "sc_expansion/csv_append.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <filesystem>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <memory>

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
         int length_cycle, std::string random_name, int random_seed, int verbosity) {

  // Staggered dimer expansion: full (non-bipartite) graph set.
  bool bipartite = false;

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    params = {Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(t_hop, 0.0), bipartite, Dual(0.0, 0.0)};
  } else {
    params = {U, beta, mu, t_hop, bipartite, 0.0};
  }

  // --- Phase 1: rank 0 generates non-bipartite vacuum diagrams, broadcasts ---
  std::vector<sc_expansion::Graph> graphs;
  if (world.rank() == 0) {
    std::cout << "Generating non-bipartite vacuum diagrams on rank 0..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    sc_expansion::VacuumDiagramGenerator gen(order, /*bipartite_only=*/false);
    gen.generate();
    graphs  = gen.get_unique_graphs();
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Generated " << graphs.size() << " unique diagrams (incl. non-bipartite) in "
              << std::chrono::duration<double>(t1 - t0).count() << " s." << std::endl;
  }
  broadcast_graphs(graphs, world);

  // --- Phase 2: dimer calculator (alpha-independent, reused across pilot + production) ---
  sc_expansion::dimer::FreeEnergyCalculator<T> calculator(params, order, graphs);

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  // --- Phase 3: alpha auto-tuning pilot ---
  // Run a short MC with the user-supplied alpha as initial guess, then rescale via
  //   alpha_new = alpha_0 * <|Omega|/W> / <alpha/W>
  // which equates alpha to typical|Omega|, the variance-minimising choice for the
  // ratio estimator coeff = alpha * beta^n * <Omega/W> / <alpha/W>.
  {
    int pilot_warmup    = 1000;
    int pilot_cycles    = 5000;
    int pilot_block     = 200;
    int pilot_n_bins    = std::max(20, pilot_cycles / pilot_block);
    int pilot_blk_size  = (pilot_cycles / pilot_n_bins) + 1;

    auto pilot_config = std::make_unique<sc_expansion::dimer::Configuration<T>>(params, order, calculator, alpha);
    triqs::mc_tools::mc_generic<double> pilot_mc(random_name, random_seed, 0);
    measure_dimer<T> pilot_meas(pilot_config.get(), pilot_n_bins, pilot_blk_size, alpha, 0);
    pilot_mc.add_move(move<T>(pilot_config.get(), pilot_mc.get_rng()), "time_swap");
    pilot_mc.add_measure(pilot_meas, "alpha_pilot");

    if (world.rank() == 0) { std::cout << "Pilot run for alpha tuning (initial alpha = " << alpha << ")..." << std::endl; }
    pilot_mc.warmup_and_accumulate(pilot_warmup, pilot_cycles, length_cycle, triqs::utility::clock_callback(-1));
    pilot_mc.collect_results(world);

    double new_alpha = alpha;
    if (world.rank() == 0) {
      double mean_abs_omega = pilot_meas.result->mean_omega_abs; // <|Omega|/W>
      double mean_alpha_W   = pilot_meas.result->mean_abs;       // <alpha/W>
      if (mean_alpha_W > 1e-18 && mean_abs_omega > 0.0) {
        new_alpha = alpha * mean_abs_omega / mean_alpha_W;
      }
      std::cout << "Pilot estimated typical|Omega| = " << new_alpha << "  (was alpha = " << alpha << ")" << std::endl;
    }
    MPI_Bcast(&new_alpha, 1, MPI_DOUBLE, 0, world.get());
    alpha = new_alpha;
  }

  // --- Phase 4: production MC with tuned alpha ---
  auto config = std::make_unique<sc_expansion::dimer::Configuration<T>>(params, order, calculator, alpha);

  triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);

  int target_block_size = 2000;
  int n_bins            = std::max(50, n_cycles / target_block_size);
  int block_size        = (n_cycles / n_bins) + 1;

  measure_dimer<T> meas(config.get(), n_bins, block_size, alpha, verbosity);
  mc.add_move(move<T>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, "dimer_coefficient");

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
    std::cout << "Order-" << order << " staggered dimer coefficient: " << meas.result->coeff << " ± " << meas.result->error << std::endl;

    std::string observable = std::is_same_v<T, Dual> ? "density" : "Omega";
    std::string filename   = "./results/dimer_square_lattice_" + observable + ".csv";

    std::ostringstream row;
    row << std::setprecision(17);
    row << U << ',' << beta << ',' << mu << ',' << t_hop << ',' << order << ',' << alpha << ',' << meas.result->coeff << ',' << meas.result->error;

    sc_expansion::append_csv_row(filename, "U,beta,mu,t,order,alpha,coeff,error", row.str());
  }
}

int main(int argc, char *argv[]) {

  if (argc < 7) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu t [alpha] [use_dual]" << std::endl;
    }
    return 1;
  }

  int order     = std::stoi(argv[1]);
  int n_cycles  = std::stoi(argv[2]);
  double U      = std::stod(argv[3]);
  double beta   = std::stod(argv[4]);
  double mu     = std::stod(argv[5]);
  double t_hop  = std::stod(argv[6]);
  double alpha  = (argc > 7 ? std::stod(argv[7]) : 0.001);
  bool use_dual = (argc > 8 ? std::stoi(argv[8]) != 0 : false);

  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    std::cout << "=== Strong Coupling MC (Staggered Dimer) ===" << std::endl;
    std::cout << "MPI ranks: " << world.size() << std::endl;
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " t=" << t_hop << " alpha=" << alpha
              << " (non-bipartite graphs)" << std::endl;
  }

  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (use_dual) {
    run<Dual>(world, order, n_cycles, U, beta, mu, t_hop, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  } else {
    run<double>(world, order, n_cycles, U, beta, mu, t_hop, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity);
  }

  return 0;
}
