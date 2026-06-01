#pragma once
#include <deque>
#include <vector>
#include "../graph.hpp"
#include "../hubbard_solver.hpp"
#include "diagram.hpp"
#include "vertex.hpp"

namespace sc_expansion::dimer {

  // Build the rooted catalog for a dimer density-density measurement at physical
  // displacement r. Runs DimerDistanceRootedDiagramGenerator and flattens its
  // (V → unique rooted graphs) output into parallel vectors of (Graph, marks,
  // sites). Each graph already carries the rooted symmetry factor and
  // free_multiplicity = 1 via the Graph override constructor; marks are in the
  // graph's canonical labeling and `sites_out[i]` gives the within-dimer site
  // (0 or 1) of each mark (the atomic catalog has no site analog). Intended as
  // the rank-0-only step of an MPI driver, with (graphs_out, marks_out,
  // sites_out) broadcast to the other ranks.
  void build_rooted_catalog(int order, std::vector<int> const &r, std::vector<Graph> &graphs_out, std::vector<std::vector<int>> &marks_out,
                            std::vector<std::vector<int>> &sites_out);

  // Diagram-list owner + evaluator for the dimer rooted density-density series.
  // Mirrors atomic::SumDiagrams but finite-U only — the dimer has no infinite-U
  // path (HilbertTraits<2> never invokes G0n_infinite_U), so all infinite_U
  // arguments, the *_infinite_U_coefficient methods, and the SJT infinite-U
  // sweep are dropped. The vacuum/free-energy path is owned separately by
  // dimer::FreeEnergyCalculator and is intentionally not duplicated here.
  //
  // Unlike atomic there is NO per-diagram scalar lattice_multiplier: on the
  // staggered superlattice the per-vertex cumulant value depends on the bond
  // directions, so the embedding multiplicity is folded into the rooted
  // Diagram's weighted spatial configurations (task 3) rather than carried as a
  // separate count. density_density(taus) is therefore just the diagram sum.
  template <typename T> class SumDiagrams {
    public:
    // Rooted density-density constructor: builds the catalog of rooted dimer
    // topologies for physical displacement r and mark spins (s1, s2), then
    // builds one rooted Diagram per catalog entry. order = truncation order.
    SumDiagrams(Parameters<T> const &params, int order, std::vector<int> r, int s1, int s2);

    // Prebuilt-catalog constructor (MPI: rank 0 enumerates via
    // build_rooted_catalog and broadcasts (graphs, marks, sites) to all ranks).
    // graphs[i] must already carry the rooted symmetry factor and
    // free_multiplicity = 1; marks[i]/sites[i] give the two marks' vertex
    // indices and within-dimer sites in graphs[i]'s labeling.
    SumDiagrams(Parameters<T> const &params, int order, std::vector<Graph> const &graphs, std::vector<std::vector<int>> const &marks,
                std::vector<std::vector<int>> const &sites, std::vector<int> const &r, int s1, int s2);

    // Finite-cluster rooted density-density constructor: builds the catalog for r
    // and mark spins (s1, s2), then builds one cluster-restricted rooted Diagram
    // per entry on the given superlattice cluster, with per-dimer normaliser
    // n_cluster_sites. Lets the same finite cluster used by the free-energy MCMC
    // carry a ⟨n(r)n(0)⟩ measurement.
    //
    // pin_origin (default false) selects the reference-site convention of the
    // cluster rooted Diagram: false = sweep mark0 over all cells ÷n_cluster_sites
    // (per-dimer average); true = pin mark0 at cluster_positions[0] with no
    // division (single-site correlator matching finite-cluster ED at that site).
    SumDiagrams(Parameters<T> const &params, int order, std::vector<int> r, int s1, int s2,
                std::vector<std::pair<int, int>> const &cluster_positions, int n_cluster_sites, bool pin_origin = false);

    // Series for the target displacement r at the given imaginary times. The
    // whole series collapses to a single scalar (the dimer targets one r vector,
    // so there is no d²-shell map). `taus` has size n_lines — the static
    // densities live at τ=0 and consume no hopping τ, so no padding slot.
    T density_density(std::vector<double> const &taus) const;

    // FreeEnergyCalculator-compatible surface, so SumDiagrams is a drop-in
    // calculator for the templated dimer::Configuration MCMC driver (whose
    // integrand f becomes the ⟨n(r)n(0)⟩ series value instead of the free
    // energy). compute_sum_diagrams is just an alias for density_density;
    // clear_all_caches drops every cached value (≡ mark_all_dirty).
    T compute_sum_diagrams(std::vector<double> const &taus) const { return this->density_density(taus); }
    void clear_all_caches() { this->mark_all_dirty(); }

    void mark_tau_dirty(int tau_index);
    void mark_all_dirty();

    std::deque<Graph> const &get_graphs() const { return this->graphs; }
    std::deque<VertexType<T>> const &get_vertex_types() const { return this->vertex_types; }
    std::deque<Diagram<T>> const &get_diagrams() const { return this->diagrams; }
    HubbardSolver<2, T> const &get_solver() const { return this->solver; }

    int get_n_diagrams() const { return (int)this->diagrams.size(); }

    bool is_density_density_mode() const { return !this->target_r.empty(); }
    std::vector<int> const &get_target_r() const { return this->target_r; }

    private:
    void init_from_rooted_catalog(std::vector<Graph> const &graphs, std::vector<std::vector<int>> const &marks,
                                  std::vector<std::vector<int>> const &sites, std::vector<int> const &r, int s1, int s2,
                                  std::vector<std::pair<int, int>> const &cluster_positions = {}, int n_cluster_sites = 0,
                                  bool pin_origin = false);

    Parameters<T> const &params;
    int order;
    std::deque<Graph> graphs;
    std::deque<Diagram<T>> diagrams;
    std::deque<VertexType<T>> vertex_types; // index = cumulant_order - 1
    HubbardSolver<2, T> solver;
    std::vector<int> target_r; // empty ⇒ not in density-density mode
  };

} // namespace sc_expansion::dimer
