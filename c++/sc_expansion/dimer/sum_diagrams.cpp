#include "sum_diagrams.hpp"
#include "../dual.hpp"
#include "../generate_diagrams.hpp"
#include <algorithm>
#include <numeric>

namespace {
  std::vector<int> sorted_graph_indices(std::vector<sc_expansion::Graph> const &graphs) {
    std::vector<int> idx(graphs.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](int a, int b) { return graphs[a].get_V() > graphs[b].get_V(); });
    return idx;
  }
} // namespace

namespace sc_expansion::dimer {

  void build_rooted_catalog(int order, std::vector<int> const &r, std::vector<Graph> &graphs_out, std::vector<std::vector<int>> &marks_out,
                            std::vector<std::vector<int>> &sites_out) {
    graphs_out.clear();
    marks_out.clear();
    sites_out.clear();

    DimerDistanceRootedDiagramGenerator gen(r, order);
    gen.generate();

    for (auto const &[V, bucket] : gen.get_rooted_graphs()) {
      (void)V;
      for (auto const &rg : bucket) {
        graphs_out.emplace_back(rg.get_canonical_form(), rg.get_V(), rg.get_rooted_automorphism_count(), (int)rg.get_rooted_symmetry_factor(),
                                /*free_multiplicity=*/1, rg.get_bipartite_only());
        marks_out.push_back(rg.get_marks());
        sites_out.push_back(rg.get_sites());
      }
    }
  }

  template <typename T>
  void SumDiagrams<T>::init_from_rooted_catalog(std::vector<Graph> const &graphs_in, std::vector<std::vector<int>> const &marks,
                                                std::vector<std::vector<int>> const &sites, std::vector<int> const &r, int s1, int s2,
                                                std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites,
                                                bool pin_origin) {
    this->target_r  = r;
    this->n_catalog = (int)graphs_in.size();
    this->n_pruned  = 0;

    // Max cumulant order over the catalog graphs, accounting for the
    // StaticDensity mark bonus (each mark on a vertex adds +1 to that vertex's
    // effective degree). This only ever over-allocates vertex_types relative to
    // what the rooted Diagram actually indexes (cumulant_order = degree/2), so
    // it is always safe; mirrors atomic's init_from_rooted_catalog.
    int max_cumulant_order = std::max(1, this->order / 2);
    for (size_t i = 0; i < graphs_in.size(); ++i) {
      auto const &g  = graphs_in[i];
      auto const &mk = marks[i];
      int V          = g.get_V();
      for (int v = 0; v < V; ++v) {
        int deg = 0;
        for (int j = 0; j < V; ++j) deg += g(v, j) + g(j, v);
        int mark_bonus = 0;
        for (int m : mk)
          if (m == v) ++mark_bonus;
        max_cumulant_order = std::max(max_cumulant_order, deg / 2 + mark_bonus);
      }
    }
    for (int k = 1; k <= max_cumulant_order; ++k) this->vertex_types.emplace_back(2 * k);

    std::vector<VertexType<T> *> vt_ptrs(max_cumulant_order);
    for (int k = 0; k < max_cumulant_order; ++k) vt_ptrs[k] = &this->vertex_types[k];

    std::vector<int> mark_spins = {s1, s2};
    auto order_idx              = sorted_graph_indices(graphs_in);
    bool on_cluster             = !cluster_positions.empty();
    for (int i : order_idx) {
      this->graphs.emplace_back(graphs_in[i]);
      if (on_cluster)
        this->diagrams.emplace_back(this->graphs.back(), vt_ptrs, marks[i], sites[i], mark_spins, r, cluster_positions, n_cluster_sites,
                                    MarkEncoding::StaticDensity, pin_origin);
      else
        this->diagrams.emplace_back(this->graphs.back(), vt_ptrs, marks[i], sites[i], mark_spins, r, MarkEncoding::StaticDensity);

      // Drop zero-contribution diagrams (mirrors atomic's embedding_counts[i]==0
      // skip in atomic::SumDiagrams::init_from_rooted_catalog). The
      // DimerDistanceRootedDiagramGenerator's d_super/diameter filters are only
      // NECESSARY conditions for ≥1 mark-pinned embedding, so the actual
      // mark-constrained enumeration (compute_spatial_configurations_rooted[_cluster])
      // can still come up empty — e.g. a graph whose only parity-allowed dimer
      // sector is unreachable at this order, or a graph that does not fit the open
      // cluster. Such a diagram has zero total embedding weight (free_multiplicity
      // == 0) and hence no valid configurations, so it contributes EXACTLY zero to
      // density_density(); keeping it only wastes per-MC-step evaluate()/dirty
      // bookkeeping. Pruning here changes no physics (the summed series is
      // identical) — it only shrinks the evaluated diagram list.
      //
      // graphs/diagrams are std::deque, so pop_back/emplace_back never invalidate
      // references to the other elements: the Graph const& that each retained
      // Diagram holds stays valid, and dropping the just-built tail pair (diagram
      // first — it references graphs.back() — then its graph) is safe.
      Diagram<T> const &built = this->diagrams.back();
      if (built.get_free_multiplicity() == 0.0 || built.get_valid_configurations().empty()) {
        this->diagrams.pop_back();
        this->graphs.pop_back();
        ++this->n_pruned;
      }
    }
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, std::vector<int> r, int s1, int s2)
     : params(params_), order(order_), solver(params_) {
    std::vector<Graph> rooted_graphs;
    std::vector<std::vector<int>> marks;
    std::vector<std::vector<int>> sites;
    build_rooted_catalog(this->order, r, rooted_graphs, marks, sites);
    this->init_from_rooted_catalog(rooted_graphs, marks, sites, r, s1, s2);
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, std::vector<Graph> const &graphs_, std::vector<std::vector<int>> const &marks,
                              std::vector<std::vector<int>> const &sites, std::vector<int> const &r, int s1, int s2)
     : params(params_), order(order_), solver(params_) {
    this->init_from_rooted_catalog(graphs_, marks, sites, r, s1, s2);
  }

  template <typename T>
  SumDiagrams<T>::SumDiagrams(Parameters<T> const &params_, int order_, std::vector<int> r, int s1, int s2,
                              std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites, bool pin_origin)
     : params(params_), order(order_), solver(params_) {
    std::vector<Graph> rooted_graphs;
    std::vector<std::vector<int>> marks;
    std::vector<std::vector<int>> sites;
    build_rooted_catalog(this->order, r, rooted_graphs, marks, sites);
    this->init_from_rooted_catalog(rooted_graphs, marks, sites, r, s1, s2, cluster_positions, n_cluster_sites, pin_origin);
  }

  template <typename T> T SumDiagrams<T>::density_density(std::vector<double> const &taus) const {
    T sum = T(0.0);
    for (auto &diagram : this->diagrams) {
      T val = const_cast<Diagram<T> &>(diagram).evaluate(taus, this->solver);
      sum   = sum + val;
    }
    return sum;
  }

  template <typename T> void SumDiagrams<T>::mark_tau_dirty(int tau_index) {
    for (auto &diagram : this->diagrams) diagram.mark_tau_dirty(tau_index);
  }

  template <typename T> void SumDiagrams<T>::mark_all_dirty() {
    for (auto &diagram : this->diagrams) diagram.mark_all_dirty();
  }

  template class SumDiagrams<double>;
  template class SumDiagrams<Dual>;

} // namespace sc_expansion::dimer
