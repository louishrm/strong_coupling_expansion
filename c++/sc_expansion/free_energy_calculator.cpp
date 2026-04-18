#include "free_energy_calculator.hpp"
#include "generate_diagrams.hpp"
#include "combinatorics.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>

#include "dual.hpp"

namespace {
  // Returns indices into `graphs` sorted by descending vertex count.
  // Diagrams with more vertices are evaluated first to seed the global
  // VertexType cache with diverse cumulants before hitting the expensive
  // concentrated diagrams (e.g. 2 vertices each with order-5 cumulants).
  std::vector<int> sorted_graph_indices(std::vector<sc_expansion::Graph> const &graphs) {
    std::vector<int> idx(graphs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return graphs[a].get_V() > graphs[b].get_V(); });
    return idx;
  }
} // namespace

namespace sc_expansion {

  template <int N_sites, typename T>
  FreeEnergyCalculator<N_sites, T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_, int override_fm_)
     : params(params_), order(order_), solver(params_) {
    VacuumDiagramGenerator gen(this->order, params.bipartite);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    // Create VertexType objects for cumulant orders 1..order/2 (max possible at any vertex)
    int max_cumulant_order = this->order / 2;
    for (int k = 1; k <= max_cumulant_order; k++) { this->vertex_types.emplace_back(2 * k); }

    std::vector<VertexType<N_sites, T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; k++) { vt_ptrs[k] = &this->vertex_types[k]; }

    // Sort by descending vertex count for optimal global cache seeding
    auto order_idx = sorted_graph_indices(unique_graphs);
    for (int i : order_idx) {
      auto const &g = unique_graphs[i];
      if (override_fm_ >= 0) {
        this->graphs.emplace_back(g.get_canonical_form(), g.get_V(), g.get_automorphism_count(), (int)g.get_symmetry_factor(), override_fm_,
                    g.get_bipartite_only());
      } else {
        this->graphs.emplace_back(g);
      }
      this->diagrams.emplace_back(this->graphs.back(), vt_ptrs);
    }
  }

  template <int N_sites, typename T>
  FreeEnergyCalculator<N_sites, T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_,
                                                         std::vector<Graph> const &prebuilt_graphs)
     : params(params_), order(order_), solver(params_) {

    int max_cumulant_order = this->order / 2;
    for (int k = 1; k <= max_cumulant_order; k++) { this->vertex_types.emplace_back(2 * k); }

    std::vector<VertexType<N_sites, T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; k++) { vt_ptrs[k] = &this->vertex_types[k]; }

    auto order_idx = sorted_graph_indices(prebuilt_graphs);
    for (int i : order_idx) {
      this->graphs.emplace_back(prebuilt_graphs[i]);
      this->diagrams.emplace_back(this->graphs.back(), vt_ptrs);
    }
  }

  template <int N_sites, typename T>
  FreeEnergyCalculator<N_sites, T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_,
                                                         std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites)
     : params(params_), order(order_), solver(params_) {
    VacuumDiagramGenerator gen(this->order, params.bipartite);
    gen.generate();
    const auto &unique_graphs = gen.get_unique_graphs();

    // Create VertexType objects for cumulant orders 1..order/2
    int max_cumulant_order = this->order / 2;
    for (int k = 1; k <= max_cumulant_order; k++) { this->vertex_types.emplace_back(2 * k); }

    std::vector<VertexType<N_sites, T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; k++) { vt_ptrs[k] = &this->vertex_types[k]; }

    auto order_idx = sorted_graph_indices(unique_graphs);
    for (int i : order_idx) {
      this->graphs.emplace_back(unique_graphs[i]);
      this->diagrams.emplace_back(this->graphs.back(), vt_ptrs, cluster_positions, n_cluster_sites);
    }
  }

  template <int N_sites, typename T>
  T FreeEnergyCalculator<N_sites, T>::compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const {
    T sum = T(0.0);
    for (auto &diagram : this->diagrams) { sum = sum + const_cast<Diagram<N_sites, T>&>(diagram).evaluate(taus, this->solver, infinite_U); }
    return sum;
  }

  template <int N_sites, typename T>
  T FreeEnergyCalculator<N_sites, T>::compute_sum_diagrams_dimer(std::vector<double> const &taus, bool infinite_U) const {
    T sum = T(0.0);
    for (auto &diagram : this->diagrams) { sum = sum + const_cast<Diagram<N_sites, T>&>(diagram).evaluate(taus, this->solver, infinite_U); }
    return sum;
  }

  template <int N_sites, typename T>
  void FreeEnergyCalculator<N_sites, T>::clear_all_caches() const {
    for (auto &vt : this->vertex_types) { vt.clear_global_cache(); }
  }

  template <int N_sites, typename T>
  void FreeEnergyCalculator<N_sites, T>::mark_tau_dirty(int tau_index) {
    for (auto &diagram : this->diagrams) { diagram.mark_tau_dirty(tau_index); }
  }

  template <int N_sites, typename T>
  void FreeEnergyCalculator<N_sites, T>::mark_all_dirty() {
    for (auto &diagram : this->diagrams) { diagram.mark_all_dirty(); }
  }

  // --- DiagMC support: single-diagram access ---

  template <int N_sites, typename T>
  T FreeEnergyCalculator<N_sites, T>::evaluate_single_diagram(int idx, std::vector<double> const &taus, bool infinite_U) const {
    return const_cast<Diagram<N_sites, T> &>(this->diagrams[idx]).evaluate(taus, this->solver, infinite_U);
  }

  template <int N_sites, typename T>
  std::vector<T> FreeEnergyCalculator<N_sites, T>::evaluate_per_config_diagram(int idx, std::vector<double> const &taus, bool infinite_U) const {
    return const_cast<Diagram<N_sites, T> &>(this->diagrams[idx]).evaluate_per_config(taus, this->solver, infinite_U);
  }

  template <int N_sites, typename T>
  void FreeEnergyCalculator<N_sites, T>::mark_single_diagram_dirty(int idx, int tau_index) {
    this->diagrams[idx].mark_tau_dirty(tau_index);
  }

  template <int N_sites, typename T>
  void FreeEnergyCalculator<N_sites, T>::mark_single_diagram_all_dirty(int idx) {
    this->diagrams[idx].mark_all_dirty();
  }

  template <int N_sites, typename T> std::pair<double, double> FreeEnergyCalculator<N_sites, T>::compute_infinite_U_coefficient(bool dimer) const {
    int n = this->order;
    std::vector<double> taus(n);
    std::iota(taus.begin(), taus.end(), 0.0);

    double sum_abs    = 0.0;
    double sum_signed = 0.0;
    do {
      // All taus change between permutations — mark every VertexInstance dirty
      for (auto &d : this->diagrams) { const_cast<Diagram<N_sites, T>&>(d).mark_all_dirty(); }
      T val_T = dimer ? this->compute_sum_diagrams_dimer(taus, true) : this->compute_sum_diagrams(taus, true);
      double val;
      if constexpr (std::is_same_v<T, Dual>) {
        val = val_T.derivative;
      } else {
        val = (double)val_T;
      }
      sum_abs += std::abs(val);
      sum_signed += val;
    } while (std::next_permutation(taus.begin(), taus.end()));

    double fact = 1.0;
    for (int i = 1; i <= n; ++i) fact *= i;

    double beta_val;
    if constexpr (std::is_same_v<T, Dual>) {
      beta_val = this->params.beta.value;
    } else {
      beta_val = (double)this->params.beta;
    }

    std::pair<double, double> result;
    result.first  = (std::pow(beta_val, n) / fact) * sum_abs;
    result.second = (std::pow(beta_val, n) / fact) * sum_signed;
    return result;
  }

  template class FreeEnergyCalculator<1, double>;
  template class FreeEnergyCalculator<1, Dual>;
  template class FreeEnergyCalculator<2, double>;
  template class FreeEnergyCalculator<2, Dual>;

} // namespace sc_expansion
