// MCMC for the staggered (triangular) dimer expansion — free energy AND the
// equal-time density-density correlator, in one file (mirrors mcmc_atom.cpp).
//
// The mode is chosen exactly as in mcmc_atom.cpp: if a displacement r=(rx,ry)
// is supplied on the command line we run the rooted ⟨n_{r,s2}(0) n_{0,s1}(0)⟩
// correlator series; otherwise we run the free-energy (Omega) series. The two
// modes share the same templated dimer::Configuration, move and defensive
// measure_dimer estimator — only the calculator changes:
//   * free energy : dimer::FreeEnergyCalculator (vacuum diagrams)
//   * correlator  : dimer::SumDiagrams (rooted density-density catalog)
//
// The dimer superlattice is non-bipartite (each dimer has 6 NN), so both
// catalogs span bipartite + non-bipartite topologies; params.bipartite is
// forced false accordingly.
//
// Both modes enumerate their diagram catalog deterministically on every rank
// (free energy: the vacuum-diagram set from `order`; correlator: the rooted
// catalog from (order, r)), so each rank builds an identical calculator
// independently — no broadcast. The vacuum enumeration is CPU-heavy at high
// order but its peak memory is only the deduped unique-graph set (~MB), so
// running it per rank is safe even with many ranks per node; the ranks would
// otherwise sit idle while rank 0 generated, so this costs no wall-clock.
//
// Correlator mode embeds on the full (infinite) lattice by default
// (use_cluster=0), or on a 3-dimer triangular cluster (use_cluster=1,
// ED-comparable).
//
// use_dual (μ-derivative → density) is supported for the free-energy expansion
// only. The density-density correlator has no μ-derivative observable, so
// use_dual is rejected in that mode (as in mcmc_atom.cpp) and the correlator
// always runs in double precision.
//
// Reference scheme (both modes): Configuration samples with weight W = |f+alpha|
// and the coefficient is recovered via the ratio estimator
//   coeff = alpha * beta^n * <f/W> / <alpha/W>,
// with f the free-energy or ⟨n(r)n(0)⟩ series value at the sampled times.

#include "sc_expansion/dimer/configuration.hpp"
#include "sc_expansion/dimer/free_energy_calculator.hpp"
#include "sc_expansion/dimer/sum_diagrams.hpp"
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

// Shared pilot + production MCMC for either calculator. `ConfigT` is the
// dimer::Configuration instantiation wrapping `Calculator` (FreeEnergyCalculator
// for the free energy, SumDiagrams for the correlator); both take the same
// (params, order, calculator, alpha) constructor. The pilot auto-tunes alpha
// in place (returned through the reference), and the measured coefficient/error
// is returned to the caller for CSV output.
template <typename T, typename ConfigT, typename Calculator>
DimerMeasureResult run_mcmc(mpi::communicator &world, sc_expansion::Parameters<T> const &params, int order, int n_cycles, Calculator &calculator,
                            double &alpha, std::string const &measure_name, int n_warmup_cycles, int length_cycle, std::string random_name,
                            int random_seed, int verbosity) {

  // --- Alpha auto-tuning pilot ---
  // Run a short MC with the user-supplied alpha as initial guess, then rescale via
  //   alpha_new = alpha_0 * <|f|/W> / <alpha/W>
  // which equates alpha to typical|f|, the variance-minimising choice for the
  // ratio estimator coeff = alpha * beta^n * <f/W> / <alpha/W>.
  {
    int pilot_warmup   = 1000;
    int pilot_cycles   = 2000;
    int pilot_block    = 200;
    int pilot_n_bins   = std::max(20, pilot_cycles / pilot_block);
    int pilot_blk_size = (pilot_cycles / pilot_n_bins) + 1;

    auto pilot_config = std::make_unique<ConfigT>(params, order, calculator, alpha);
    triqs::mc_tools::mc_generic<double> pilot_mc(random_name, random_seed, 0);
    measure_dimer<T> pilot_meas(pilot_config.get(), pilot_n_bins, pilot_blk_size, alpha, 0);
    pilot_mc.add_move(move<T>(pilot_config.get(), pilot_mc.get_rng()), "time_swap");
    pilot_mc.add_measure(pilot_meas, "alpha_pilot");

    if (world.rank() == 0) { std::cout << "Pilot run for alpha tuning (initial alpha = " << alpha << ")..." << std::endl; }
    pilot_mc.warmup_and_accumulate(pilot_warmup, pilot_cycles, length_cycle, triqs::utility::clock_callback(-1));
    pilot_mc.collect_results(world);

    double new_alpha = alpha;
    if (world.rank() == 0) {
      double mean_abs_f   = pilot_meas.result->mean_omega_abs; // <|f|/W>
      double mean_alpha_W = pilot_meas.result->mean_abs;       // <alpha/W>
      if (mean_alpha_W > 1e-18 && mean_abs_f > 0.0) { new_alpha = alpha * mean_abs_f / mean_alpha_W; }
      std::cout << "Pilot estimated typical|f| = " << new_alpha << "  (was alpha = " << alpha << ")" << std::endl;
    }
    MPI_Bcast(&new_alpha, 1, MPI_DOUBLE, 0, world.get());
    alpha = new_alpha;
  }

  // --- Production MC with tuned alpha ---
  auto config = std::make_unique<ConfigT>(params, order, calculator, alpha);

  triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);

  int target_block_size = 2000;
  int n_bins            = std::max(50, n_cycles / target_block_size);
  int block_size        = (n_cycles / n_bins) + 1;

  measure_dimer<T> meas(config.get(), n_bins, block_size, alpha, verbosity);
  mc.add_move(move<T>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(meas, measure_name);

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

  return *meas.result;
}

// Free-energy (ln Z) expansion. Templated on the scalar T so the Dual variant
// yields a derivative coefficient: dual_mode 1 seeds μ (density ∂Ω/∂μ), dual_mode
// 2 seeds U (double occupancy ⟨n↑n↓⟩ = ∂Ω/∂U). dual_mode 0 is plain free energy.
template <typename T>
void run_free_energy(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, double t_hop, double alpha,
                     int n_warmup_cycles, int length_cycle, std::string random_name, int random_seed, int verbosity, int dual_mode) {

  // Staggered dimer expansion: full (non-bipartite) topology set.
  bool bipartite = false;

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    double mu_seed = (dual_mode == 1) ? 1.0 : 0.0;
    double U_seed  = (dual_mode == 2) ? 1.0 : 0.0;
    params         = {Dual(U, U_seed), Dual(beta, 0.0), Dual(mu, mu_seed), Dual(t_hop, 0.0), bipartite, Dual(0.0, 0.0)};
  } else {
    params = {U, beta, mu, t_hop, bipartite, 0.0};
  }

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  // Vacuum-diagram catalog. Every rank obtains an identical set with no MPI
  // exchange: it first tries the on-disk cache (warmed by generate_dimer_
  // diagrams, or self-warmed below), and on a miss falls back to generating
  // locally. The enumeration is deterministic in `order`, so the per-rank
  // fallback yields the same catalog on every rank; its peak memory is just the
  // deduped unique-graph set (~MB), so it is safe even at high rank density.
  //
  // This replaces the former rank-0-only generate-then-broadcast, which stalled
  // every other rank in a collective for the whole (~minutes at high order)
  // enumeration and then fired thousands of tiny per-graph broadcasts at once —
  // a pattern that triggered UCX wireup / connection-refused failures at large
  // rank counts.
  std::vector<sc_expansion::Graph> graphs;
  if (sc_expansion::load_vacuum_graphs(order, /*bipartite_only=*/false, graphs)) {
    if (world.rank() == 0) {
      std::cout << "Loaded " << graphs.size() << " vacuum diagrams from cache " << sc_expansion::vacuum_diagrams_path(order, false) << "."
                << std::endl;
    }
  } else {
    if (world.rank() == 0) { std::cout << "Cache miss; generating non-bipartite vacuum diagrams per rank..." << std::endl; }
    auto t0 = std::chrono::high_resolution_clock::now();
    sc_expansion::VacuumDiagramGenerator gen(order, /*bipartite_only=*/false);
    gen.generate();
    graphs  = gen.get_unique_graphs();
    auto t1 = std::chrono::high_resolution_clock::now();
    if (world.rank() == 0) {
      std::cout << "Generated " << graphs.size() << " unique diagrams (incl. non-bipartite) in " << std::chrono::duration<double>(t1 - t0).count()
                << " s." << std::endl;
      // Self-warm the cache so subsequent jobs / μ-points / dual_modes skip the
      // enumeration. Only rank 0 writes; the atomic temp+rename in save makes
      // concurrent writers (several array tasks racing a cold cache) safe.
      try {
        sc_expansion::save_vacuum_graphs(order, /*bipartite_only=*/false, graphs);
        std::cout << "Cached vacuum diagrams to " << sc_expansion::vacuum_diagrams_path(order, false) << "." << std::endl;
      } catch (std::exception const &e) { std::cerr << "Warning: could not write vacuum-diagram cache: " << e.what() << std::endl; }
    }
  }

  sc_expansion::dimer::FreeEnergyCalculator<T> calculator(params, order, graphs);

  using ConfigT = sc_expansion::dimer::Configuration<T>;
  auto result   = run_mcmc<T, ConfigT>(world, params, order, n_cycles, calculator, alpha, "dimer_coefficient", n_warmup_cycles, length_cycle,
                                       random_name, random_seed, verbosity);

  if (world.rank() == 0) {
    std::cout << "Order-" << order << " staggered dimer coefficient: " << result.coeff << " ± " << result.error << std::endl;

    std::string observable = (dual_mode == 1) ? "density" : (dual_mode == 2) ? "double_occupancy" : "free_energy";
    std::string filename   = "./results/dimer_square_lattice_" + observable + ".csv";

    std::ostringstream row;
    row << std::setprecision(17);
    row << U << ',' << beta << ',' << mu << ',' << t_hop << ',' << order << ',' << alpha << ',' << result.coeff << ',' << result.error;
    sc_expansion::append_csv_row(filename, "U,beta,mu,t,order,alpha,coeff,error", row.str());
  }
}

// Equal-time density-density correlator ⟨n(r)n(0)⟩. Double precision only:
// there is no μ-derivative observable here, so use_dual is rejected upstream.
void run_correlator(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, double t_hop, double alpha, bool use_cluster,
                    int n_warmup_cycles, int length_cycle, std::string random_name, int random_seed, int verbosity, std::vector<int> const &r, int s1,
                    int s2) {

  // Staggered dimer expansion: full (non-bipartite) topology set.
  bool bipartite                          = false;
  sc_expansion::Parameters<double> params = {U, beta, mu, t_hop, bipartite, 0.0};

  // 3-dimer triangular cluster on the staggered superlattice. n_cluster_sites is
  // the per-dimer normaliser for the sweep convention; the correlator pins one
  // reference site (pin_origin=true) so it is unused there (no ÷n_cluster_sites).
  std::vector<std::pair<int, int>> cluster_positions = {{0, 0}, {1, 0}, {0, 1}};
  int n_cluster_sites                                = 3;

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  // --- Rooted density-density calculator ---
  // Built per-rank: the catalog is a deterministic function of (order, r), so
  // every rank enumerates the same SumDiagrams independently (no broadcast).
  // Cluster mode pins n(0) at cluster_positions[0] (pin_origin=true): the
  // correlator is INTENSIVE/local, anchored at one reference site rather than
  // swept-and-averaged like the EXTENSIVE free energy, so the coefficient is
  // directly comparable to single-site finite-cluster ED. The infinite-lattice
  // mode already pins the origin in its full-lattice embedding.
  std::unique_ptr<sc_expansion::dimer::SumDiagrams<double>> calculator;
  if (use_cluster) {
    calculator = std::make_unique<sc_expansion::dimer::SumDiagrams<double>>(params, order, r, s1, s2, cluster_positions, n_cluster_sites,
                                                                            /*pin_origin=*/true);
  } else {
    calculator = std::make_unique<sc_expansion::dimer::SumDiagrams<double>>(params, order, r, s1, s2);
  }

  if (world.rank() == 0) {
    std::cout << "Building rooted catalog for r=(" << r[0] << "," << r[1] << ") spins=(" << s1 << "," << s2 << ")"
              << (use_cluster ? " on the 3-dimer cluster..." : " on the infinite lattice...") << std::endl;
    std::cout << "Rooted catalog: " << calculator->get_n_catalog() << " topologies generated, " << calculator->get_n_pruned()
              << " dropped (empty embedding), " << calculator->get_n_diagrams() << " diagrams instantiated." << std::endl;
  }

  using ConfigT = sc_expansion::dimer::Configuration<double, sc_expansion::dimer::SumDiagrams<double>>;
  auto result   = run_mcmc<double, ConfigT>(world, params, order, n_cycles, *calculator, alpha, "dimer_correlator_coefficient", n_warmup_cycles,
                                            length_cycle, random_name, random_seed, verbosity);

  if (world.rank() == 0) {
    std::cout << "Order-" << order << " dimer ⟨n(r)n(0)⟩ coefficient (r=(" << r[0] << "," << r[1] << "), spins=(" << s1 << "," << s2
              << ")): " << result.coeff << " ± " << result.error << std::endl;

    std::string filename = "./results/dimer_correlator_nn.csv";
    std::ostringstream row;
    row << std::setprecision(17);
    row << U << ',' << beta << ',' << mu << ',' << t_hop << ',' << order << ',' << r[0] << ',' << r[1] << ',' << s1 << ',' << s2 << ','
        << (use_cluster ? 1 : 0) << ',' << alpha << ',' << result.coeff << ',' << result.error;
    sc_expansion::append_csv_row(filename, "U,beta,mu,t,order,rx,ry,s1,s2,cluster,alpha,coeff,error", row.str());
  }
}

int main(int argc, char *argv[]) {

  if (argc < 7) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu t [alpha] [dual_mode] [use_cluster] [rx ry s1 s2]" << std::endl;
      std::cerr << "  rx ry s1 s2 : if all four are present, run the density-density correlator at r=(rx,ry) with mark spins (s1,s2)." << std::endl;
      std::cerr << "  dual_mode (default 0): 0 = free energy Omega, 1 = density (∂Ω/∂μ), 2 = double occupancy (∂Ω/∂U). Free-energy mode only." << std::endl;
      std::cerr << "  use_cluster (default 0): 1 = 3-dimer triangle (ED-comparable), 0 = infinite lattice (correlator mode only)." << std::endl;
      std::cerr << "  s1 s2 : mark spins (0=down, 1=up) at (0,0) and r." << std::endl;
    }
    return 1;
  }

  int order        = std::stoi(argv[1]);
  int n_cycles     = std::stoi(argv[2]);
  double U         = std::stod(argv[3]);
  double beta      = std::stod(argv[4]);
  double mu        = std::stod(argv[5]);
  double t_hop     = std::stod(argv[6]);
  double alpha     = (argc > 7 ? std::stod(argv[7]) : 0.001);
  int dual_mode    = (argc > 8 ? std::stoi(argv[8]) : 0);
  bool use_cluster = (argc > 9 ? std::stoi(argv[9]) != 0 : false);

  if (dual_mode < 0 || dual_mode > 2) {
    if (mpi::communicator().rank() == 0) { std::cerr << "dual_mode must be 0 (free energy), 1 (density), or 2 (double occupancy)." << std::endl; }
    return 1;
  }

  std::vector<int> r;
  int s1 = 0, s2 = 0;
  if (argc > 10) {
    if (argc < 14) {
      if (mpi::communicator().rank() == 0) {
        std::cerr << "Density-density mode requires all four of rx ry s1 s2; got " << (argc - 10) << "." << std::endl;
      }
      return 1;
    }
    int rx = std::stoi(argv[10]);
    int ry = std::stoi(argv[11]);
    s1     = std::stoi(argv[12]);
    s2     = std::stoi(argv[13]);
    r      = {rx, ry};
    if (dual_mode != 0) {
      if (mpi::communicator().rank() == 0) {
        std::cerr << "dual_mode must be 0 in density-density mode (no μ/U-derivative for the correlator)." << std::endl;
      }
      return 1;
    }
  }

  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    std::cout << "=== Strong Coupling MC (Staggered Dimer) ===" << std::endl;
    std::cout << "MPI ranks: " << world.size() << std::endl;
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " t=" << t_hop << " alpha=" << alpha;
    if (!r.empty())
      std::cout << " r=(" << r[0] << "," << r[1] << ") s1=" << s1 << " s2=" << s2 << (use_cluster ? " [3-dimer cluster]" : " [infinite lattice]");
    std::cout << " (non-bipartite graphs)" << std::endl;
  }

  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (!r.empty()) {
    run_correlator(world, order, n_cycles, U, beta, mu, t_hop, alpha, use_cluster, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity,
                   r, s1, s2);
  } else if (dual_mode != 0) {
    run_free_energy<Dual>(world, order, n_cycles, U, beta, mu, t_hop, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity,
                          dual_mode);
  } else {
    run_free_energy<double>(world, order, n_cycles, U, beta, mu, t_hop, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity,
                            dual_mode);
  }

  return 0;
}
