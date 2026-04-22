#include "free_energy_calculator.hpp"
#include "generate_diagrams.hpp"
#include "combinatorics.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>

#include "dual.hpp"

namespace {
  // Returns indices into `graphs` sorted by descending vertex count.
  std::vector<int> sorted_graph_indices(std::vector<sc_expansion::Graph> const &graphs) {
    std::vector<int> idx(graphs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return graphs[a].get_V() > graphs[b].get_V(); });
    return idx;
  }
} // namespace

namespace sc_expansion {

  template <typename T>
  void FreeEnergyCalculator<T>::init_from_graphs(std::vector<Graph> const &source_graphs, int override_fm) {
    // Infer max cumulant order from actual vertex degrees — works for both self-loop and pure-hopping graphs.
    int max_cumulant_order = this->order / 2;
    for (auto const &g : source_graphs) {
      int V = g.get_V();
      for (int v = 0; v < V; v++) {
        int deg = 0;
        for (int j = 0; j < V; j++) { deg += g(v, j) + g(j, v); }
        max_cumulant_order = std::max(max_cumulant_order, deg / 2);
      }
    }
    for (int k = 1; k <= max_cumulant_order; k++) { this->vertex_types.emplace_back(2 * k); }

    std::vector<VertexType<T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; k++) { vt_ptrs[k] = &this->vertex_types[k]; }

    auto order_idx = sorted_graph_indices(source_graphs);
    for (int i : order_idx) {
      auto const &g = source_graphs[i];
      if (override_fm >= 0) {
        this->graphs.emplace_back(g.get_canonical_form(), g.get_V(), g.get_automorphism_count(), (int)g.get_symmetry_factor(), override_fm,
                                  g.get_bipartite_only());
      } else {
        this->graphs.emplace_back(g);
      }
      this->diagrams.emplace_back(this->graphs.back(), vt_ptrs);
    }

    this->precompute_self_loop_diagram();
  }

  template <typename T> void FreeEnergyCalculator<T>::precompute_self_loop_diagram() {
    std::vector<double> dummy_taus(this->order, 0.0);
    for (auto &d : this->diagrams) {
      auto const &g = d.get_graph();
      if (g.get_V() == 1 && g.get_n_self_loops() == this->order) {
        this->self_loop_val_finite   = d.evaluate(dummy_taus, this->solver, false);
        this->self_loop_val_infinite = d.evaluate(dummy_taus, this->solver, true);
        this->self_loop_present      = true;
        break;
      }
    }
  }

  template <typename T> T FreeEnergyCalculator<T>::get_self_loop_contribution_finite() const {
    T contrib = this->self_loop_val_finite;
    for (int i = 0; i < this->order; i++) contrib = contrib * this->params.delta;
    for (int i = 0; i < this->order; i++) contrib = contrib * this->params.beta;
    return contrib;
  }

  template <typename T> T FreeEnergyCalculator<T>::get_self_loop_contribution_infinite() const {
    T contrib = this->self_loop_val_infinite;
    for (int i = 0; i < this->order; i++) contrib = contrib * this->params.delta;
    for (int i = 0; i < this->order; i++) contrib = contrib * this->params.beta;
    return contrib;
  }

  template <typename T>
  FreeEnergyCalculator<T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_, int override_fm_, bool allow_self_loops_)
     : params(params_), order(order_), solver(params_) {
    VacuumDiagramGenerator gen(this->order, params.bipartite, allow_self_loops_);
    gen.generate();
    this->init_from_graphs(gen.get_unique_graphs(), override_fm_);
  }

  template <typename T>
  FreeEnergyCalculator<T>::FreeEnergyCalculator(Parameters<T> const &params_, int order_, std::vector<Graph> const &prebuilt_graphs, int override_fm_)
     : params(params_), order(order_), solver(params_) {
    this->init_from_graphs(prebuilt_graphs, override_fm_);
  }

  template <typename T> T FreeEnergyCalculator<T>::compute_sum_diagrams(std::vector<double> const &taus, bool infinite_U) const {
    T sum = T(0.0);
    for (auto &diagram : this->diagrams) {
      int k = diagram.get_graph().get_n_self_loops();

      T val;
      if (this->self_loop_present && k == this->order && diagram.get_graph().get_V() == 1) {
        // Use precomputed constant — cumulant of number operators is tau-independent.
        val = infinite_U ? this->self_loop_val_infinite : this->self_loop_val_finite;
      } else {
        val = const_cast<Diagram<T> &>(diagram).evaluate(taus, this->solver, infinite_U);
      }

      T factor = T((this->order - k) % 2 == 0 ? 1.0 : -1.0);
      if (k > 0) {
        for (int i = 0; i < k; i++) { factor = factor * this->params.delta; }
      }
      sum = sum + val * factor;
    }
    return sum;
  }

  template <typename T> void FreeEnergyCalculator<T>::mark_tau_dirty(int tau_index) {
    for (auto &diagram : this->diagrams) { diagram.mark_tau_dirty(tau_index); }
  }

  template <typename T> void FreeEnergyCalculator<T>::mark_all_dirty() {
    for (auto &diagram : this->diagrams) { diagram.mark_all_dirty(); }
  }

  template <typename T> std::pair<double, double> FreeEnergyCalculator<T>::compute_infinite_U_coefficient() {
    int n = this->order;
    std::vector<double> taus(n);

    SJT sjt(n);
    double sum_abs    = 0.0;
    double sum_signed = 0.0;

    do {
      const auto &perm = sjt.get_permutation();
      for (int j = 0; j < n; j++) { taus[j] = (double)(perm[j] - 1); }

      this->mark_all_dirty();
      T val_T = this->compute_sum_diagrams(taus, true);

      double val;
      if constexpr (std::is_same_v<T, Dual>) {
        val = val_T.derivative;
      } else {
        val = (double)val_T;
      }
      sum_abs += std::abs(val);
      sum_signed += val;
    } while (sjt.next_permutation());

    double beta_val;
    if constexpr (std::is_same_v<T, Dual>) {
      beta_val = this->params.beta.value;
    } else {
      beta_val = (double)this->params.beta;
    }

    double fact = 1.0;
    for (int i = 1; i <= n; ++i) fact *= i;
    double prefactor = std::pow(beta_val, n) / fact;

    return {prefactor * sum_abs, prefactor * sum_signed};
  }

  template class FreeEnergyCalculator<double>;
  template class FreeEnergyCalculator<Dual>;

} // namespace sc_expansion
