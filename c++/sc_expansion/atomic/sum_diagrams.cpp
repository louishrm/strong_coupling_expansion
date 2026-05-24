#include "sum_diagrams.hpp"
#include "../generate_diagrams.hpp"
#include "../combinatorics.hpp"
#include "../dual.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>

namespace {
  std::vector<int> sorted_graph_indices(std::vector<sc_expansion::Graph> const &graphs) {
    std::vector<int> idx(graphs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return graphs[a].get_V() > graphs[b].get_V(); });
    return idx;
  }
} // namespace

namespace sc_expansion::atomic {

  template <typename T> void SumDiagrams<T>::init_from_graphs(std::vector<Graph> const &source_graphs, int override_fm) {
    int max_cumulant_order = this->order / 2;
    for (auto const &g : source_graphs) {
      int V = g.get_V();
      for (int v = 0; v < V; ++v) {
        int deg = 0;
        for (int j = 0; j < V; ++j) deg += g(v, j) + g(j, v);
        max_cumulant_order = std::max(max_cumulant_order, deg / 2);
      }
    }
    for (int k = 1; k <= max_cumulant_order; ++k) this->vertex_types.emplace_back(2 * k);

    std::vector<VertexType<T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; ++k) vt_ptrs[k] = &this->vertex_types[k];

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
  }

  void build_rooted_catalog(int order, bool bipartite, std::vector<int> const &r,
                            std::vector<Graph> &graphs_out,
                            std::vector<std::vector<int>> &marks_out) {
    graphs_out.clear();
    marks_out.clear();

    DistanceRootedDiagramGenerator gen(r, order, /*bipartite_only=*/bipartite);
    gen.generate();

    for (auto const &[V, bucket] : gen.get_rooted_graphs()) {
      (void)V;
      for (auto const &rg : bucket) {
        graphs_out.emplace_back(rg.get_canonical_form(), rg.get_V(), rg.get_rooted_automorphism_count(),
                                (int)rg.get_rooted_symmetry_factor(), /*free_multiplicity=*/1, rg.get_bipartite_only());
        marks_out.push_back(rg.get_marks());
      }
    }
  }

  template <typename T>
  void SumDiagrams<T>::init_from_rooted_catalog(std::vector<Graph> const &rooted_graphs, std::vector<std::vector<int>> const &marks,
                                                std::vector<int> const &r, int s1, int s2, int override_lm) {
    int d_sq = 0;
    for (int c : r) d_sq += c * c;
    this->target_d_sq = d_sq;

    // Per-diagram lattice multiplier: count of Z² embeddings of (graph, marks)
    // with mark[0] at origin and mark[1] at r. Default path computes this via
    // count_lattice_embeddings() — the DistanceRootedDiagramGenerator filter
    // is only a necessary condition for ≥1 embedding, so some catalog entries
    // have count == 0 and we drop them outright (they would contribute zero
    // anyway, just at the cost of wasted Diagram evaluation in the MC).
    //
    // override_lm ≥ 0 bypasses this: every catalog entry is kept with the
    // override as its multiplier, mirroring the "per-diagram, no lattice sum"
    // convention used by cluster-cell reference computations.
    std::vector<int> embedding_counts(rooted_graphs.size(), 0);
    if (override_lm >= 0) {
      std::fill(embedding_counts.begin(), embedding_counts.end(), override_lm);
    } else {
      for (size_t i = 0; i < rooted_graphs.size(); ++i) {
        embedding_counts[i] = count_lattice_embeddings(rooted_graphs[i], marks[i], r);
      }
    }

    // Max cumulant order over kept graphs, accounting for mark bonus (each
    // mark contributes +1 to its vertex degree under the StaticDensity
    // encoding).
    int max_cumulant_order = std::max(1, this->order / 2);
    for (size_t i = 0; i < rooted_graphs.size(); ++i) {
      if (embedding_counts[i] == 0) continue;
      auto const &g  = rooted_graphs[i];
      auto const &mk = marks[i];
      int V          = g.get_V();
      for (int v = 0; v < V; ++v) {
        int deg = 0;
        for (int j = 0; j < V; ++j) deg += g(v, j) + g(j, v);
        int mark_bonus = 0;
        for (int m : mk)
          if (m == v) ++mark_bonus;
        int co             = deg / 2 + mark_bonus;
        max_cumulant_order = std::max(max_cumulant_order, co);
      }
    }
    for (int k = 1; k <= max_cumulant_order; ++k) this->vertex_types.emplace_back(2 * k);

    std::vector<VertexType<T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; ++k) vt_ptrs[k] = &this->vertex_types[k];

    auto order_idx              = sorted_graph_indices(rooted_graphs);
    std::vector<int> mark_spins = {s1, s2};
    for (int i : order_idx) {
      if (embedding_counts[i] == 0) continue;
      this->graphs.emplace_back(rooted_graphs[i]);
      this->diagrams.emplace_back(this->graphs.back(), vt_ptrs, marks[i], mark_spins,
                                  /*flip_mark_order=*/false, MarkEncoding::StaticDensity);
      this->lattice_multiplier.push_back({{d_sq, embedding_counts[i]}});
    }
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, int override_fm_) : params(params_), order(order_), solver(params_) {
    VacuumDiagramGenerator gen(this->order, params.bipartite);
    gen.generate();
    this->init_from_graphs(gen.get_unique_graphs(), override_fm_);
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, std::vector<Graph> const &prebuilt_graphs, int override_fm_)
     : params(params_), order(order_), solver(params_) {
    this->init_from_graphs(prebuilt_graphs, override_fm_);
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, std::vector<int> r, int s1, int s2, int override_lm)
     : params(params_), order(order_), solver(params_) {
    std::vector<Graph> rooted_graphs;
    std::vector<std::vector<int>> marks;
    build_rooted_catalog(this->order, this->params.bipartite, r, rooted_graphs, marks);
    this->init_from_rooted_catalog(rooted_graphs, marks, r, s1, s2, override_lm);
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, std::vector<Graph> const &prebuilt_rooted_graphs,
                              std::vector<std::vector<int>> const &prebuilt_marks, std::vector<int> const &r, int s1, int s2, int override_lm)
     : params(params_), order(order_), solver(params_) {
    this->init_from_rooted_catalog(prebuilt_rooted_graphs, prebuilt_marks, r, s1, s2, override_lm);
  }

  template <typename T> T SumDiagrams<T>::free_energy(std::vector<double> const &taus, bool infinite_U) const {
    T sum = T(0.0);
    for (auto &diagram : this->diagrams) {
      T val = const_cast<Diagram<T> &>(diagram).evaluate(taus, this->solver, infinite_U);
      sum   = sum + val;
    }
    return sum;
  }

  template <typename T> std::map<int, T> SumDiagrams<T>::density_density(std::vector<double> const &taus, bool infinite_U) const {
    // Rooted Diagram::evaluate expects taus of size n_lines+1 with the trailing
    // entry pinned to 0 (mark-slot placeholder under StaticDensity encoding).
    std::vector<double> taus_padded;
    taus_padded.reserve(taus.size() + 1);
    taus_padded.insert(taus_padded.end(), taus.begin(), taus.end());
    taus_padded.push_back(0.0);

    std::map<int, T> sums;
    for (size_t i = 0; i < this->diagrams.size(); ++i) {
      T temporal = const_cast<Diagram<T> &>(this->diagrams[i]).evaluate(taus_padded, this->solver, infinite_U);
      for (auto const &[d_sq, mult] : this->lattice_multiplier[i]) {
        T contribution = T((double)mult) * temporal;
        auto it        = sums.find(d_sq);
        if (it == sums.end()) sums.emplace(d_sq, contribution);
        else
          it->second = it->second + contribution;
      }
    }
    return sums;
  }

  template <typename T> void SumDiagrams<T>::mark_tau_dirty(int tau_index) {
    for (auto &diagram : this->diagrams) diagram.mark_tau_dirty(tau_index);
  }

  template <typename T> void SumDiagrams<T>::mark_all_dirty() {
    for (auto &diagram : this->diagrams) diagram.mark_all_dirty();
  }

  template <typename T>
  template <class Accum, class PerPerm>
  void SumDiagrams<T>::sjt_sweep(Accum &accum, PerPerm &&per_perm) {
    int n = this->order;
    std::vector<double> taus(n);
    SJT sjt(n);
    do {
      auto const &perm = sjt.get_permutation();
      for (int j = 0; j < n; ++j) taus[j] = (double)(perm[j] - 1);
      this->mark_all_dirty();
      per_perm(taus, accum);
    } while (sjt.next_permutation());
  }

  template <typename T> std::pair<double, double> SumDiagrams<T>::free_energy_infinite_U_coefficient() {
    struct {
      double abs    = 0.0;
      double signed_ = 0.0;
    } accum;

    this->sjt_sweep(accum, [this](std::vector<double> const &taus, auto &a) {
      T val_T = this->free_energy(taus, true);
      double val;
      if constexpr (std::is_same_v<T, Dual>) {
        val = val_T.derivative;
      } else {
        val = (double)val_T;
      }
      a.abs += std::abs(val);
      a.signed_ += val;
    });

    double beta_val;
    if constexpr (std::is_same_v<T, Dual>) {
      beta_val = this->params.beta.value;
    } else {
      beta_val = (double)this->params.beta;
    }
    double fact = 1.0;
    for (int i = 1; i <= this->order; ++i) fact *= i;
    double prefactor = std::pow(beta_val, this->order) / fact;
    return {prefactor * accum.abs, prefactor * accum.signed_};
  }

  template <typename T> std::map<int, std::pair<double, double>> SumDiagrams<T>::density_density_infinite_U_coefficient() {
    std::map<int, std::pair<double, double>> accum;

    this->sjt_sweep(accum, [this](std::vector<double> const &taus, auto &a) {
      auto sums = this->density_density(taus, true);
      for (auto const &[d_sq, val_T] : sums) {
        double val;
        if constexpr (std::is_same_v<T, Dual>) {
          val = val_T.derivative;
        } else {
          val = (double)val_T;
        }
        auto &pr = a[d_sq];
        pr.first += std::abs(val);
        pr.second += val;
      }
    });

    double beta_val;
    if constexpr (std::is_same_v<T, Dual>) {
      beta_val = this->params.beta.value;
    } else {
      beta_val = (double)this->params.beta;
    }
    double fact = 1.0;
    for (int i = 1; i <= this->order; ++i) fact *= i;
    double prefactor = std::pow(beta_val, this->order) / fact;

    for (auto &kv : accum) {
      kv.second.first *= prefactor;
      kv.second.second *= prefactor;
    }
    return accum;
  }

  template class SumDiagrams<double>;
  template class SumDiagrams<Dual>;

} // namespace sc_expansion::atomic
