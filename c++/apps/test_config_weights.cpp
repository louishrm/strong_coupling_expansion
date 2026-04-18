// Diagnostic: probe how the per-diagram integrand decomposes across its
// valid global configurations, for tau vectors drawn from the Metropolis
// chain used by the full-sum dimer MC (weight = |Ω - Ω_∞U| + α |Ω_∞U|).
//
// Pipeline:
//   1. Build DimerConfiguration + triqs mc_generic, run n_warmup cycles.
//   2. Accumulate: every `length_cycle` MC steps, snapshot config.state.
//   3. For each snapshot tau, compute per-config signed contributions for
//      every diagram and aggregate statistics.
//
// Metrics per diagram (averaged over snapshots):
//   - cancellation  = <|sum c_i|> / <sum |c_i|>
//   - N_eff         = <(sum |c_i|)^2 / sum c_i^2>
//   - max/L1        = <max |c_i| / sum |c_i|>
//   - topK%         = L1 fraction captured by top K% of configs (sorted per sample)
//
// Usage:
//   mpirun -n <ranks> test_config_weights <order> <n_samples> <U> <beta> <mu> <t_hop> [alpha=1.0] [n_warmup=1000] [length_cycle=10]

#include "sc_expansion/dimer_configuration.hpp"
#include "sc_expansion/move.hpp"
#include "sc_expansion/diagram.hpp"
#include "sc_expansion/generate_diagrams.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>
#include <mpi/mpi.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace sc_expansion;

static void broadcast_graphs(std::vector<Graph> &graphs, mpi::communicator &world) {
  int n_graphs = (int)graphs.size();
  MPI_Bcast(&n_graphs, 1, MPI_INT, 0, world.get());
  if (world.rank() != 0) { graphs.clear(); }
  for (int g = 0; g < n_graphs; g++) {
    int V, aut, sym, fm, bip_only;
    std::vector<uint8_t> adj;
    if (world.rank() == 0) {
      auto const &gr = graphs[g];
      V              = gr.get_V();
      aut            = gr.get_automorphism_count();
      sym            = (int)gr.get_symmetry_factor();
      fm             = (int)gr.get_free_multiplicity();
      bip_only       = gr.get_bipartite_only() ? 1 : 0;
      adj            = gr.get_canonical_form();
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

// Measure that just snapshots tau state each cycle.
struct tau_snapshot_measure {
  DimerConfiguration<double> *config;
  std::vector<std::vector<double>> *snapshots;

  tau_snapshot_measure(DimerConfiguration<double> *c, std::vector<std::vector<double>> *s) : config(c), snapshots(s) {}

  void accumulate(double /*sign*/) { this->snapshots->push_back(this->config->state); }
  void collect_results(mpi::communicator const & /*c*/) {}
};

struct DiagramStats {
  long n_samples        = 0;
  int n_configs         = 0;
  int V                 = 0;
  double sum_L1         = 0.0;
  double sum_abs_signed = 0.0;
  double sum_neff       = 0.0;
  double sum_max_frac   = 0.0;
  std::vector<double> sum_topk_frac;
};

static const std::vector<double> PERCENTILES = {0.01, 0.05, 0.10, 0.25, 0.50, 0.90};

int main(int argc, char *argv[]) {
  if (argc < 7) {
    if (mpi::communicator().rank() == 0) {
      std::cerr << "Usage: " << argv[0] << " order n_samples U beta mu t_hop [alpha=1.0] [n_warmup=1000] [length_cycle=10]" << std::endl;
    }
    return 1;
  }
  int order        = std::stoi(argv[1]);
  long n_samples   = std::stol(argv[2]);
  double U         = std::stod(argv[3]);
  double beta      = std::stod(argv[4]);
  double mu        = std::stod(argv[5]);
  double t_hop     = std::stod(argv[6]);
  double alpha     = (argc > 7 ? std::stod(argv[7]) : 1.0);
  int n_warmup     = (argc > 8 ? std::stoi(argv[8]) : 100);
  int length_cycle = (argc > 9 ? std::stoi(argv[9]) : 1);

  mpi::environment env(argc, argv);
  mpi::communicator world;
  int rank = world.rank(), size = world.size();

  long n_local = n_samples / size + (rank < (n_samples % size) ? 1 : 0);

  if (rank == 0) {
    std::cout << "=== Config-weight diagnostic (dimer, MC-drawn taus) ===" << std::endl;
    std::cout << "order=" << order << " U=" << U << " beta=" << beta << " mu=" << mu << " t=" << t_hop << " alpha=" << alpha << " ranks=" << size
              << " samples_total=" << n_samples << " warmup=" << n_warmup << " length_cycle=" << length_cycle << std::endl;
  }

  std::vector<Graph> graphs;
  if (rank == 0) {
    VacuumDiagramGenerator gen(order, true);
    gen.generate();
    graphs = gen.get_unique_graphs();
    std::cout << "Generated " << graphs.size() << " unique diagrams" << std::endl;
  }
  broadcast_graphs(graphs, world);

  Parameters<double> params{U, beta, mu, t_hop, true};
  FreeEnergyCalculator<2, double> calculator(params, order, graphs);

  auto const &diagrams = calculator.get_diagrams();
  int n_diag           = (int)diagrams.size();

  // --- MCMC warmup + sampling ---
  auto config = std::make_unique<DimerConfiguration<double>>(params, order, calculator, alpha);
  triqs::mc_tools::mc_generic<double> mc("", 32186222 + rank * 786512, rank == 0 ? 2 : 0);

  std::vector<std::vector<double>> snapshots;
  snapshots.reserve((size_t)n_local);

  mc.add_move(move<double>(config.get(), mc.get_rng()), "time_swap");
  mc.add_measure(tau_snapshot_measure(config.get(), &snapshots), "snapshot");

  auto t_mc0 = std::chrono::high_resolution_clock::now();
  mc.warmup_and_accumulate((int)n_warmup, (int)n_local, length_cycle, triqs::utility::clock_callback(-1));
  auto t_mc1 = std::chrono::high_resolution_clock::now();

  if (rank == 0) {
    std::cout << "MC done in " << std::chrono::duration<double>(t_mc1 - t_mc0).count() << " s (rank 0); " << snapshots.size() << " snapshots/rank"
              << std::endl;
  }

  // --- Per-snapshot diagnostic ---
  std::vector<DiagramStats> stats(n_diag);
  for (int d = 0; d < n_diag; d++) {
    stats[d].n_configs = (int)diagrams[d].get_valid_configurations().size();
    stats[d].V         = diagrams[d].get_graph().get_V();
    stats[d].sum_topk_frac.assign(PERCENTILES.size(), 0.0);
  }

  auto t0 = std::chrono::high_resolution_clock::now();
  for (auto const &taus : snapshots) {
    for (int d = 0; d < n_diag; d++) {
      auto contribs = calculator.evaluate_per_config_diagram(d, taus, false);

      double L1 = 0.0, signed_sum = 0.0, L2sq = 0.0, maxabs = 0.0;
      std::vector<double> abs_contribs;
      abs_contribs.reserve(contribs.size());
      for (auto c : contribs) {
        double a = std::abs(c);
        abs_contribs.push_back(a);
        L1 += a;
        L2sq += c * c;
        signed_sum += c;
        if (a > maxabs) maxabs = a;
      }
      if (L1 == 0.0) {
        stats[d].n_samples++;
        continue;
      }

      std::sort(abs_contribs.begin(), abs_contribs.end(), std::greater<double>());
      double cum = 0.0;
      int n      = (int)abs_contribs.size();
      int pi     = 0;
      for (int k = 0; k < n && pi < (int)PERCENTILES.size(); k++) {
        cum += abs_contribs[k];
        while (pi < (int)PERCENTILES.size() && (double)(k + 1) / n >= PERCENTILES[pi]) {
          stats[d].sum_topk_frac[pi] += cum / L1;
          pi++;
        }
      }
      while (pi < (int)PERCENTILES.size()) {
        stats[d].sum_topk_frac[pi] += 1.0;
        pi++;
      }

      stats[d].sum_L1 += L1;
      stats[d].sum_abs_signed += std::abs(signed_sum);
      stats[d].sum_neff += (L1 * L1) / L2sq;
      stats[d].sum_max_frac += maxabs / L1;
      stats[d].n_samples++;
    }
  }
  auto t1           = std::chrono::high_resolution_clock::now();
  double analysis_s = std::chrono::duration<double>(t1 - t0).count();

  std::vector<DiagramStats> global = stats;
  for (int d = 0; d < n_diag; d++) {
    long ns_g   = 0;
    double L1_g = 0, abs_g = 0, neff_g = 0, maxf_g = 0;
    MPI_Reduce(&stats[d].n_samples, &ns_g, 1, MPI_LONG, MPI_SUM, 0, world.get());
    MPI_Reduce(&stats[d].sum_L1, &L1_g, 1, MPI_DOUBLE, MPI_SUM, 0, world.get());
    MPI_Reduce(&stats[d].sum_abs_signed, &abs_g, 1, MPI_DOUBLE, MPI_SUM, 0, world.get());
    MPI_Reduce(&stats[d].sum_neff, &neff_g, 1, MPI_DOUBLE, MPI_SUM, 0, world.get());
    MPI_Reduce(&stats[d].sum_max_frac, &maxf_g, 1, MPI_DOUBLE, MPI_SUM, 0, world.get());
    std::vector<double> tk_g(PERCENTILES.size(), 0.0);
    MPI_Reduce(stats[d].sum_topk_frac.data(), tk_g.data(), (int)PERCENTILES.size(), MPI_DOUBLE, MPI_SUM, 0, world.get());
    if (rank == 0) {
      global[d].n_samples      = ns_g;
      global[d].sum_L1         = L1_g;
      global[d].sum_abs_signed = abs_g;
      global[d].sum_neff       = neff_g;
      global[d].sum_max_frac   = maxf_g;
      global[d].sum_topk_frac  = tk_g;
    }
  }

  if (rank == 0) {
    std::cout << "\nAnalysis time (s, rank 0): " << analysis_s << std::endl;
    std::cout << "\n--- Per-diagram config distribution (averaged over " << n_samples << " MC-drawn tau samples) ---" << std::endl;
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "d   V Ncfg    cancel   N_eff  max/L1 |";
    for (double p : PERCENTILES) std::cout << "  top" << std::setprecision(0) << std::setw(3) << (p * 100) << "%";
    std::cout << std::setprecision(4) << std::endl;
    std::cout << std::string(40 + 9 * (int)PERCENTILES.size(), '-') << std::endl;
    for (int d = 0; d < n_diag; d++) {
      auto const &g = global[d];
      double ns     = (double)g.n_samples;
      double cancel = g.sum_abs_signed / g.sum_L1;
      double neff   = g.sum_neff / ns;
      double maxf   = g.sum_max_frac / ns;
      std::cout << std::setw(3) << d << " " << std::setw(2) << g.V << " " << std::setw(5) << g.n_configs << "  " << std::setw(7) << cancel << "  "
                << std::setw(7) << std::setprecision(2) << neff << "  " << std::setprecision(4) << std::setw(6) << maxf << " |";
      for (size_t i = 0; i < PERCENTILES.size(); i++) { std::cout << "  " << std::setw(6) << g.sum_topk_frac[i] / ns; }
      std::cout << std::endl;
    }
    std::cout << "\nReading the table:" << std::endl;
    std::cout << "  cancel   = <|sum c|>/<sum|c|>  (1 = no cancellation)" << std::endl;
    std::cout << "  N_eff    = effective # configs contributing (smaller = more peaked)" << std::endl;
    std::cout << "  max/L1   = fraction held by single largest config" << std::endl;
    std::cout << "  topK%    = L1 fraction captured by the top K% of configs (sorted per sample)" << std::endl;
  }
  return 0;
}
