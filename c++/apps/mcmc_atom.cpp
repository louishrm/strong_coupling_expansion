#include "sc_expansion/atomic/configuration.hpp"
#include "sc_expansion/atomic/sum_diagrams.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/combinatorics.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/dimer/measure_dimer.hpp"
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
#include <limits>
#include <cstdlib>

// Broadcast a rooted catalog (graphs + per-graph mark pair) from rank 0 to all
// other ranks. Each graph is serialised as broadcast_graphs() does; the marks
// pair travels alongside as two extra ints per graph. The receiving rank
// reconstructs (graphs, marks) using Graph's override constructor — i.e. the
// caller is responsible for having passed rooted symmetry factor / fm = 1 in.
static void broadcast_rooted_catalog(std::vector<sc_expansion::Graph> &graphs, std::vector<std::vector<int>> &marks, mpi::communicator &world) {

  int n_graphs = (int)graphs.size();
  MPI_Bcast(&n_graphs, 1, MPI_INT, 0, world.get());

  if (world.rank() != 0) {
    graphs.clear();
    marks.clear();
  }

  for (int g = 0; g < n_graphs; g++) {
    int V, aut, sym, fm;
    int bip_only;
    int mark0, mark1;
    std::vector<uint8_t> adj;

    if (world.rank() == 0) {
      auto const &graph = graphs[g];
      V                 = graph.get_V();
      aut               = graph.get_automorphism_count();
      sym               = (int)graph.get_symmetry_factor();
      fm                = (int)graph.get_free_multiplicity();
      bip_only          = graph.get_bipartite_only() ? 1 : 0;
      adj               = graph.get_canonical_form();
      mark0             = marks[g][0];
      mark1             = marks[g][1];
    }

    MPI_Bcast(&V, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&aut, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&sym, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&fm, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&bip_only, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&mark0, 1, MPI_INT, 0, world.get());
    MPI_Bcast(&mark1, 1, MPI_INT, 0, world.get());

    int adj_size = V * V;
    if (world.rank() != 0) { adj.resize(adj_size); }
    MPI_Bcast(adj.data(), adj_size, MPI_UNSIGNED_CHAR, 0, world.get());

    if (world.rank() != 0) {
      graphs.emplace_back(adj, V, aut, sym, fm, bip_only != 0);
      marks.push_back({mark0, mark1});
    }
  }
}

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
//
// Since the atomic MC moved to the uniform-reference scheme (W = |f + alpha|),
// this is computed only for diagnostic reporting (printed + written to CSV) and
// does not feed the estimator.
template <typename T>
std::pair<double, double> compute_reference_integral_mpi(sc_expansion::atomic::SumDiagrams<T> &calculator,
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

  bool corr_mode = calculator.is_density_density_mode();

  for (uint64_t i = 0; i < my_count; i++) {
    auto const &perm = sjt.get_permutation();
    for (int j = 0; j < order; j++) { taus[j] = (double)(perm[j] - 1); }

    calculator.mark_all_dirty();
    T val_T;
    if (corr_mode) {
      val_T = calculator.density_density_single(taus, true);
    } else {
      val_T = calculator.free_energy(taus, true);
    }
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
         int length_cycle, std::string random_name, int random_seed, int verbosity, std::vector<int> const &r, int s1, int s2, int dual_mode,
         bool compute_reference) {

  bool corr_mode = !r.empty();

  sc_expansion::Parameters<T> params;
  if constexpr (std::is_same_v<T, Dual>) {
    // dual_mode picks which static parameter carries the unit derivative seed:
    //   1 → μ  (∂Ω/∂μ = density)
    //   2 → U  (∂Ω/∂U = double occupancy ⟨n↑n↓⟩, by Hellmann-Feynman on the U·n↑n↓ term)
    double mu_seed = (dual_mode == 1) ? 1.0 : 0.0;
    double U_seed  = (dual_mode == 2) ? 1.0 : 0.0;
    params         = {Dual(U, U_seed), Dual(beta, 0.0), Dual(mu, mu_seed), Dual(0.0, 0.0), bipartite, Dual(0.0, 0.0)};
  } else {
    params = {U, beta, mu, 0.0, bipartite, 0.0};
  }

  // --- Phase 1: Rank 0 builds the catalog, then broadcasts to all ranks ---
  std::vector<sc_expansion::Graph> graphs;
  std::vector<std::vector<int>> rooted_marks;
  if (world.rank() == 0) {
    if (corr_mode) {
      std::cout << "Generating rooted (density-density) catalog on rank 0 for r=(";
      for (size_t i = 0; i < r.size(); ++i) std::cout << r[i] << (i + 1 < r.size() ? "," : "");
      std::cout << ")..." << std::endl;
      auto t0 = std::chrono::high_resolution_clock::now();
      sc_expansion::atomic::build_rooted_catalog(order, params.bipartite, r, graphs, rooted_marks);
      auto t1 = std::chrono::high_resolution_clock::now();
      std::cout << "Generated " << graphs.size() << " rooted topologies in " << std::chrono::duration<double>(t1 - t0).count() << " s." << std::endl;
    } else {
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
  }
  if (corr_mode) {
    broadcast_rooted_catalog(graphs, rooted_marks, world);
  } else {
    broadcast_graphs(graphs, world);
  }

  // --- Phase 2: All ranks construct the calculator from pre-built graphs ---
  std::unique_ptr<sc_expansion::atomic::SumDiagrams<T>> calculator_ptr;
  if (corr_mode) {
    calculator_ptr = std::make_unique<sc_expansion::atomic::SumDiagrams<T>>(params, order, graphs, rooted_marks, r, s1, s2);
  } else {
    calculator_ptr = std::make_unique<sc_expansion::atomic::SumDiagrams<T>>(params, order, graphs);
  }
  auto &calculator = *calculator_ptr;

  // --- Phase 3: MPI-parallel infinite-U reference integral using SJT (diagnostic only) ---
  // Bypassed by default (compute_reference=false). The production estimator is the
  // uniform-reference |f+alpha| scheme (Phases 4-5), which needs no infinite-U
  // baseline; the SJT integral is an O(n!) diagnostic only — and the double-occupancy
  // (dual-U) observable has no infinite-U baseline to compare against in the first
  // place. Set SCE_COMPUTE_REFERENCE=1 to re-enable it (e.g. to validate the Omega
  // series against the exact infinite-U coefficient).
  double reference_integral        = std::numeric_limits<double>::quiet_NaN();
  double signed_reference_integral = std::numeric_limits<double>::quiet_NaN();
  if (compute_reference) {
    if (world.rank() == 0) { std::cout << "Computing reference integral across " << world.size() << " MPI ranks (SJT) [diagnostic only]..." << std::endl; }
    auto t_ref_start = std::chrono::high_resolution_clock::now();
    auto [ri, sri]   = compute_reference_integral_mpi(calculator, params, order, world);
    auto t_ref_end   = std::chrono::high_resolution_clock::now();
    reference_integral        = ri;
    signed_reference_integral = sri;
    if (world.rank() == 0) {
      std::cout << "Reference integral: " << signed_reference_integral << " (abs: " << reference_integral << ")"
                << " computed in " << std::chrono::duration<double>(t_ref_end - t_ref_start).count() << " s." << std::endl;
    }
  } else if (world.rank() == 0) {
    std::cout << "Skipping infinite-U reference integral (|f+alpha| estimator only; set SCE_COMPUTE_REFERENCE=1 to enable the diagnostic)." << std::endl;
  }

  if (world.rank() == 0) { std::filesystem::create_directory("./results"); }

  // --- Phase 4: Alpha auto-tuning pilot ---
  // Short MC with user-supplied alpha as initial guess; rescale via
  //   alpha_new = alpha_0 * <|f|/W> / <alpha/W>
  // so alpha matches typical|f|, the variance-minimising choice for the ratio
  // estimator coeff = alpha * beta^n * <f/W> / <alpha/W>.
  {
    int pilot_warmup   = 1000;
    int pilot_cycles   = 5000;
    int pilot_block    = 200;
    int pilot_n_bins   = std::max(20, pilot_cycles / pilot_block);
    int pilot_blk_size = (pilot_cycles / pilot_n_bins) + 1;

    auto pilot_config = std::make_unique<sc_expansion::atomic::Configuration<T>>(params, order, alpha, calculator);
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

  // --- Phase 5: Production MC with tuned alpha ---
  auto config = std::make_unique<sc_expansion::atomic::Configuration<T>>(params, order, alpha, calculator);

  triqs::mc_tools::mc_generic<double> mc(random_name, random_seed, verbosity);

  int target_block_size = 2000;
  int n_bins            = std::max(50, n_cycles / target_block_size);
  int block_size        = (n_cycles / n_bins) + 1;

  measure_dimer<T> meas(config.get(), n_bins, block_size, alpha, verbosity);
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
    if (compute_reference) { std::cout << "Exact infinite-U coefficient (order " << order << "): " << signed_reference_integral << std::endl; }

    std::string lattice = bipartite ? "square" : "triangular";

    namespace fs = std::filesystem;
    std::string filename;
    std::string header;
    std::string row_str;
    if (corr_mode) {
      filename = "./results/atom_" + lattice + "_lattice_corr.csv";
      header   = "U,beta,mu,order,alpha,rx,ry,s1,s2,coeff,error,reference_integral";
      std::ostringstream row;
      row << std::setprecision(17);
      row << U << ',' << beta << ',' << mu << ',' << order << ',' << alpha << ',' << r[0] << ',' << r[1] << ',' << s1 << ',' << s2 << ','
          << meas.result->coeff << ',' << meas.result->error << ',' << reference_integral;
      row_str = row.str();
    } else {
      std::string observable = (dual_mode == 1) ? "density" : (dual_mode == 2) ? "double_occupancy" : "free_energy";
      filename               = "./results/atom_" + lattice + "_lattice_" + observable + ".csv";
      header                 = "U,beta,mu,order,alpha,coeff,error,reference_integral";
      std::ostringstream row;
      row << std::setprecision(17);
      row << U << ',' << beta << ',' << mu << ',' << order << ',' << alpha << ',' << meas.result->coeff << ',' << meas.result->error << ','
          << reference_integral;
      row_str = row.str();
    }

    std::cout << "Writing row to " << fs::absolute(filename).string() << " ..." << std::endl;
    sc_expansion::append_csv_row(filename, header, row_str);
    std::cout << "CSV write OK (file size " << fs::file_size(filename) << " bytes)." << std::endl;
  }
}

int main(int argc, char *argv[]) {

  if (argc < 6) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_cycles U beta mu [bipartite] [alpha] [dual_mode] [rx ry s1 s2]" << std::endl;
      std::cerr << "  dual_mode (default 0): 0 = free energy (Omega); 1 = density (∂Ω/∂μ); 2 = double occupancy (∂Ω/∂U)." << std::endl;
      std::cerr << "  rx ry s1 s2: if all four are present, run density-density correlator at r=(rx,ry) with mark spins (s1,s2)." << std::endl;
      std::cerr << "  env SCE_COMPUTE_REFERENCE=1: also compute the (diagnostic, O(n!)) infinite-U reference integral; default skips it." << std::endl;
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
  int dual_mode  = (argc > 8 ? std::stoi(argv[8]) : 0); // 0=Omega, 1=density (∂/∂μ), 2=double occupancy (∂/∂U)
  bool use_dual  = (dual_mode != 0);

  if (dual_mode < 0 || dual_mode > 2) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "dual_mode must be 0 (Omega), 1 (density, ∂/∂μ) or 2 (double occupancy, ∂/∂U); got " << dual_mode << "." << std::endl;
    }
    return 1;
  }

  // Diagnostic infinite-U reference integral is bypassed by default (the |f+alpha|
  // estimator is self-contained); set SCE_COMPUTE_REFERENCE=1 to re-enable it.
  const char *ref_env    = std::getenv("SCE_COMPUTE_REFERENCE");
  bool compute_reference = (ref_env != nullptr && std::string(ref_env) == "1");

  std::vector<int> r;
  int s1 = 0, s2 = 0;
  if (argc > 9) {
    if (argc < 13) {
      if (mpi::communicator().rank() == 0) {
        std::cerr << "Density-density mode requires all four of rx ry s1 s2; got " << (argc - 9) << "." << std::endl;
      }
      return 1;
    }
    int rx = std::stoi(argv[9]);
    int ry = std::stoi(argv[10]);
    s1     = std::stoi(argv[11]);
    s2     = std::stoi(argv[12]);
    r      = {rx, ry};
    if (use_dual) {
      if (mpi::communicator().rank() == 0) {
        std::cerr << "dual_mode must be 0 in density-density mode (no μ/U-derivative for the correlator)." << std::endl;
      }
      return 1;
    }
  }

  mpi::environment env(argc, argv);
  mpi::communicator world;

  if (world.rank() == 0) {
    std::cout << "=== Strong Coupling MC (Atomic) ===" << std::endl;
    std::cout << "MPI ranks: " << world.size() << std::endl;
    std::cout << "Order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " bipartite=" << bipartite << " alpha=" << alpha
              << " dual_mode=" << dual_mode;
    if (!r.empty()) std::cout << " r=(" << r[0] << "," << r[1] << ") s1=" << s1 << " s2=" << s2;
    std::cout << std::endl;
  }

  int length_cycle        = 1;
  int n_warmup_cycles     = 2000;
  std::string random_name = "";
  int random_seed         = 32186222 + world.rank() * 786512;
  int verbosity           = (world.rank() == 0 ? 2 : 0);

  if (use_dual) {
    run<Dual>(world, order, n_cycles, U, beta, mu, bipartite, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity, r, s1, s2,
              dual_mode, compute_reference);
  } else {
    run<double>(world, order, n_cycles, U, beta, mu, bipartite, alpha, n_warmup_cycles, length_cycle, random_name, random_seed, verbosity, r, s1, s2,
                dual_mode, compute_reference);
  }

  return 0;
}
