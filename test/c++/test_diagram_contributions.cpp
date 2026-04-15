/*
 * Per-diagram contribution analysis on the infinite lattice (N_sites=2 dimer).
 *
 * For orders 4 and 6, evaluates each diagram independently at random tau
 * points and prints the fraction of the total integrand carried by each
 * diagram, grouped by vertex count.
 *
 * This answers: how important is the V=2 diagram on the infinite lattice
 * (where high-V diagrams have large free multiplicities)?
 *
 * Usage:  ./test_diagram_contributions   (single rank, no MPI needed)
 */

#include <gtest/gtest.h>
#include "sc_expansion/diagram.hpp"
#include "sc_expansion/graph.hpp"
#include "sc_expansion/free_energy_calculator.hpp"
#include <random>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <map>
#include <set>

using namespace sc_expansion;

static void run_contribution_analysis(int order, Parameters<double> const &params, int n_samples) {

  FreeEnergyCalculator<2, double> calculator(params, order);
  HubbardSolver<2, double> solver(params);

  auto const &diags  = calculator.get_diagrams();
  auto const &graphs = calculator.get_graphs();
  int n_diags        = (int)diags.size();

  double beta = params.beta;
  std::mt19937 rng(12345 + order);
  std::uniform_real_distribution<double> tau_dist(0.0, beta);

  std::vector<double> avg_abs(n_diags, 0.0);
  std::vector<double> time_per_diag(n_diags, 0.0);

  for (int s = 0; s < n_samples; s++) {
    std::vector<double> taus(order);
    for (int i = 0; i < order; i++) taus[i] = tau_dist(rng);

    for (int d = 0; d < n_diags; d++) {
      auto &diagram = const_cast<Diagram<2, double> &>(diags[d]);
      diagram.mark_all_dirty();
      auto t0    = std::chrono::high_resolution_clock::now();
      double val = diagram.evaluate(taus, solver, false);
      auto t1    = std::chrono::high_resolution_clock::now();
      avg_abs[d] += std::abs(val);
      time_per_diag[d] += std::chrono::duration<double>(t1 - t0).count();
    }
  }

  // Normalize
  double sum_all = 0, sum_time = 0;
  for (int d = 0; d < n_diags; d++) {
    avg_abs[d] /= n_samples;
    sum_all += avg_abs[d];
    sum_time += time_per_diag[d];
  }

  // Print per-diagram
  std::cout << "\n=== Per-diagram contributions (Order " << order << ", infinite lattice, " << n_samples << " samples) ===" << std::endl;
  std::cout << std::left << std::setw(6) << "Diag" << std::setw(6) << "V" << std::setw(10) << "SymF" << std::setw(18) << "<|D_d|>"
            << std::setw(12) << "fraction" << std::setw(12) << "time (s)" << std::setw(12) << "time %" << std::endl;
  std::cout << std::string(76, '-') << std::endl;

  for (int d = 0; d < n_diags; d++) {
    double frac      = 100.0 * avg_abs[d] / sum_all;
    double time_frac = 100.0 * time_per_diag[d] / sum_time;
    std::cout << std::left << std::setw(6) << d << std::setw(6) << graphs[d].get_V() << std::setw(10) << graphs[d].get_symmetry_factor()
              << std::setw(18) << std::scientific << std::setprecision(6) << avg_abs[d] << std::fixed << std::setprecision(2) << std::setw(12) << frac
              << std::setw(12) << time_per_diag[d] << std::setw(12) << time_frac << "%" << std::endl;
  }

  // Group by V
  std::map<int, double> by_V, time_by_V;
  for (int d = 0; d < n_diags; d++) {
    by_V[graphs[d].get_V()] += avg_abs[d];
    time_by_V[graphs[d].get_V()] += time_per_diag[d];
  }

  std::cout << "\n--- Grouped by vertex count ---" << std::endl;
  std::cout << std::left << std::setw(6) << "V" << std::setw(10) << "n_diags" << std::setw(18) << "sum <|D_d|>" << std::setw(12) << "signal %"
            << std::setw(12) << "time (s)" << std::setw(12) << "time %" << std::endl;
  std::cout << std::string(70, '-') << std::endl;

  std::map<int, int> count_by_V;
  for (int d = 0; d < n_diags; d++) count_by_V[graphs[d].get_V()]++;

  for (auto &[V, val] : by_V) {
    double frac      = 100.0 * val / sum_all;
    double time_frac = 100.0 * time_by_V[V] / sum_time;
    std::cout << std::left << std::setw(6) << V << std::setw(10) << count_by_V[V] << std::setw(18) << std::scientific << std::setprecision(6) << val
              << std::fixed << std::setprecision(2) << std::setw(12) << frac << std::setw(12) << time_by_V[V] << std::setw(12) << time_frac << "%"
              << std::endl;
  }
  std::cout << "\nTotal sum <|D_d|> = " << std::scientific << std::setprecision(6) << sum_all << std::endl;
  std::cout << "Total time = " << std::fixed << std::setprecision(2) << sum_time << " s" << std::endl;
  std::cout << std::string(60, '=') << "\n" << std::endl;
}

// =====================================================================
// Diagnostic: count distinct local states under vertex-level symmetries
//
// For each diagram, apply spin flip (op ^ 2) and site swap (op ^ 1)
// to each local_state and check how many orbits remain.
// =====================================================================

static std::vector<uint8_t> apply_to_local_state(std::vector<uint8_t> const &ops, int symmetry) {
  // symmetry: 0 = identity, 1 = spin flip, 2 = site swap, 3 = both
  std::vector<uint8_t> result(ops.size());
  for (size_t i = 0; i < ops.size(); i++) {
    uint8_t o = ops[i];
    if (symmetry & 1) o ^= 2; // spin flip: flip bit 1
    if (symmetry & 2) o ^= 1; // site swap: flip bit 0
    result[i] = o;
  }
  return result;
}

static int count_orbits(std::vector<std::vector<uint8_t>> const &states, int symmetry_mask) {
  // symmetry_mask: bitmask of which symmetries to apply
  //   bit 0 = spin flip (op ^ 2), bit 1 = site swap (op ^ 1)
  // Returns the number of distinct orbits under the symmetry group.
  std::set<std::vector<uint8_t>> canonical;
  for (auto const &ops : states) {
    // Generate the orbit of this state
    std::vector<uint8_t> min_rep = ops;
    for (int sym = 0; sym < 4; sym++) {
      if ((sym & ~symmetry_mask) != 0) continue; // skip symmetries not in mask
      auto transformed = apply_to_local_state(ops, sym);
      if (transformed < min_rep) min_rep = transformed;
    }
    canonical.insert(min_rep);
  }
  return (int)canonical.size();
}

static void run_local_state_symmetry_analysis(int order, Parameters<double> const &params) {

  FreeEnergyCalculator<2, double> calculator(params, order);
  auto const &diags  = calculator.get_diagrams();
  auto const &graphs = calculator.get_graphs();

  std::cout << "\n=== Local state symmetry analysis (Order " << order << ", infinite lattice) ===" << std::endl;
  std::cout << std::left << std::setw(6) << "Diag" << std::setw(6) << "V" << std::setw(10) << "SymF"
            << std::setw(10) << "Vertex" << std::setw(14) << "RawStates" << std::setw(14) << "SpinFlip"
            << std::setw(14) << "SiteSwap" << std::setw(14) << "Both" << std::endl;
  std::cout << std::string(88, '-') << std::endl;

  for (size_t d = 0; d < diags.size(); d++) {
    auto const &all_ls = diags[d].get_local_states();
    int V = graphs[d].get_V();

    for (int v = 0; v < V; v++) {
      auto const &states = all_ls[v];
      int raw       = (int)states.size();
      int spin_only = count_orbits(states, 1);  // spin flip only
      int site_only = count_orbits(states, 2);  // site swap only
      int both      = count_orbits(states, 3);  // both

      // Only print vertices with significant state counts
      if (raw < 10) continue;

      std::cout << std::left << std::setw(6) << d << std::setw(6) << V << std::setw(10) << graphs[d].get_symmetry_factor()
                << std::setw(10) << v << std::setw(14) << raw << std::setw(14) << spin_only
                << std::setw(14) << site_only << std::setw(14) << both << std::endl;
    }
  }
  std::cout << std::string(88, '=') << "\n" << std::endl;
}

TEST(DiagramContributions, Order4InfiniteLattice) {
  Parameters<double> params{8.0, 2.0, 3.0, 1.0, true};
  run_contribution_analysis(4, params, 5000);
}

TEST(DiagramContributions, Order6InfiniteLattice) {
  Parameters<double> params{8.0, 2.0, 3.0, 1.0, true};
  run_contribution_analysis(6, params, 2000);
}

TEST(DiagramContributions, Order4LocalStateSymmetry) {
  Parameters<double> params{8.0, 2.0, 3.0, 1.0, true};
  run_local_state_symmetry_analysis(4, params);
}

TEST(DiagramContributions, Order6LocalStateSymmetry) {
  Parameters<double> params{8.0, 2.0, 3.0, 1.0, true};
  run_local_state_symmetry_analysis(6, params);
}

TEST(DiagramContributions, Order8LocalStateSymmetry) {
  Parameters<double> params{8.0, 2.0, 3.0, 1.0, true};
  run_local_state_symmetry_analysis(8, params);
}
