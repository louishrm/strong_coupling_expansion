#include "sc_expansion/atomic/configuration.hpp"
#include "sc_expansion/atomic/free_energy_calculator.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/combinatorics.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/measure.hpp"
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

// Broadcast a vector of Graph objects from rank 0 to all other ranks.
// Each graph is serialised as: V, automorphism_count, symmetry_factor, free_multiplicity,
// bipartite_only flag, and the flattened canonical adjacency matrix (V*V uint8_t values).
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

// MPI-parallel computation of the infinite-U reference integral using SJT.
// The n! permutations (simplices) are divided into contiguous chunks in SJT order.
// Each rank fast-forwards the SJT generator to its chunk start, then evaluates its
// portion. Within each chunk the consecutive-transposition property is preserved,
// so the vertex cache is kept alive across permutations for maximum reuse.
template <typename T>
std::pair<double, double> compute_reference_integral_mpi(sc_expansion::atomic::FreeEnergyCalculator<T> &calculator,
                                                         sc_expansion::Parameters<T> const &params, int order, mpi::communicator &world) {

  uint64_t n_perms = sc_expansion::factorial(order);
  uint64_t rank    = world.rank();
  uint64_t size    = world.size();

  uint64_t chunk_size = n_perms / size;
  uint64_t remainder  = n_perms % size;

  uint64_t my_start = rank * chunk_size + std::min(rank, remainder);
  uint64_t my_count = chunk_size + (rank < remainder ? 1 : 0);

  sc_expansion::SJT sjt(order);

  for (uint64_t i = 0; i < my_start; i++) { sjt.next_permutation(); }

  double local_sum_abs    = 0.0;
  double local_sum_signed = 0.0;

  std::vector<double> taus(order);

  for (uint64_t i = 0; i < my_count; i++) {
    auto const &perm = sjt.get_permutation();
    for (int j = 0; j < order; j++) { taus[j] = (double)(perm[j] - 1); }

    calculator.mark_all_dirty();
    T val_T = calculator.compute_sum_diagrams(taus, true);
    double val;
    if constexpr (std::is_same_v<T, Dual>) {
      val = val_T.derivative;
    } else {
      val = (double)val_T;
    }
    local_sum_abs += std::abs(val);
    local_sum_signed += val;

    if (i + 1 < my_count) { sjt.next_permutation(); }
  }

  double global_sum_abs    = 0.0;
  double global_sum_signed = 0.0;
  MPI_Allreduce(&local_sum_abs, &global_sum_abs, 1, MPI_DOUBLE, MPI_SUM, world.get());
  MPI_Allreduce(&local_sum_signed, &global_sum_signed, 1, MPI_DOUBLE, MPI_SUM, world.get());

  double beta_val;
  if constexpr (std::is_same_v<T, Dual>) {
    beta_val = params.beta.value;
  } else {
    beta_val = (double)params.beta;
  }
  double fact = 1.0;
  for (int i = 1; i <= order; ++i) fact *= i;

  double norm = std::pow(beta_val, order) / fact;
  return {norm * global_sum_abs, norm * global_sum_signed};
}

template <typename T>
void run(mpi::communicator &world, int order, int n_cycles, double U, double beta, double mu, bool bipartite, double alpha, int n_warmup_cycles,
         int length_cycle, std::string random_name, int random_seed, int verbosity) {

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    params = {Dual(U, 0.0), Dual(beta, 0.0), Dual(mu, 1.0), Dual(0.0, 0.0), bipartite, Dual(0.0, 0.0)};
  } else {
    params = {U, beta, mu, 0.0, bipartite, 0.0};
  }

  // --- Phase 1: Rank 0 generates all vacuum diagrams, then broadcasts to all ranks ---
  std::vector<sc_expansion::Graph> graphs;
  if (world.rank() == 0) {
    bool loaded = false;
    if (params.bipartite) {
      auto path = sc_expansion::bipartite_diagrams_path(order);
      if (sc_expansion::load_bipartite_graphs(order, graphs)) {
        std::cout << "Loaded " << graphs.size() << " cached diagrams from " << path << std::endl;
        loaded = true;
      }
    }
    if (!loaded) {
      std::cout << "Generating vacuum diagrams on rank 0..." << std::endl;
      auto t0 = std::chrono::high_resolution_clock::now();
      sc_expansion::VacuumDiagramGenerator gen(order, params.bipartite);
      gen.generate();
      graphs  = gen.get_unique_graphs();
      auto t1 = std::chrono::high_resolution_clock::now();
      std::cout << "Generated " << graphs.size() << " unique diagrams in " << std::chrono::duration<double>(t1 - t0).count() << " s." << std::endl;
    }
  }
  broadcast_graphs(graphs, world);

  // --- Phase 2: All ranks construct the calculator from pre-built graphs ---
  sc_expansion::atomic::FreeEnergyCalculator<T> calculator(params, order, graphs);

  // --- Phase 3: MPI-parallel infinite-U reference integral using SJT ---
  if (world.rank() == 0) { std::cout << "Computing reference integral across " << world.size() << " MPI ranks (SJT)..." << std::endl; }
  auto t_ref_start                                     = std::chrono::high_resolution_clock::now();
  auto [reference_integral, signed_reference_integral] = compute_reference_integral_mpi(calculator, params, order, world);
  auto t_ref_end                                       = std::chrono::high_resolution_clock::now();
  if (world.rank() == 0) {
    std::cout << "Reference integral: " << signed_reference_integral << " (abs: " << reference_integral << ")"
              << " computed in " << std::chrono::duration<double>(t_ref_end - t_ref_start).count() << " s." << std::endl;
  }

  // --- Phase 4: MC sampling (reuses the same calculator — no second diagram generation) ---
  auto config = std::make_unique<sc_expansion::atomic::Configuration<T>>(params, order, alpha, calculator);

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);

  int target_block_size = 2000;
  int n_bins            = std::max(50, n_cycles / target_block_size);
  int block_size        = (n_cycles / n_bins) + 1;

  measure<T> meas(config.get(), reference_integral, signed_reference_integral, n_bins, block_size, mu, verbosity);
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

    std::string lattice    = bipartite ? "square" : "triangular";
    std::string observable = std::is_same_v<T, Dual> ? "density" : "Omega";
    std::string filename   = "./results/atom_" + lattice + "_lattice_" + observable + ".csv";

    std::ostringstream row;
    row << std::setprecision(17);
    row << U << ',' << beta << ',' << mu << ',' << order << ',' << alpha << ',' << meas.result->mean << ',' << meas.result->error << ','
        << reference_integral;

    sc_expansion::append_csv_row(filename, "U,beta,mu,order,alpha,coeff,error,reference_integral", row.str());
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
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " bipartite=" << bipartite << " alpha=" << alpha
              << std::endl;
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
