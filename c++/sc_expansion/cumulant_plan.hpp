#pragma once
#include <vector>
#include <cstdint>
#include "args.hpp"
#include "hubbard_solver.hpp"

namespace sc_expansion {

  // Precomputed τ-independent decomposition plan for one vertex.
  // Node i holds:
  //   - a leaf G0 call (by global indices into the vertex's master op list)
  //   - a list of product-terms to ADD to the G0 leaf (signs already include all combinatorial signs)
  // Evaluation is a single forward sweep:
  //   value[i] = G0^(k)(leaf) + Σ term.sign * Π value[term.factor_node_ids]
  // The final answer is value[root_id].
  //
  // The plan is plain data (no T, no N_sites) — only integer indices and signs.
  struct CumulantPlan {

    struct LeafOps {
      // Ordered global indices (into the master unprimed / primed op arrays of the vertex).
      // Order matches the bitmask ctz-iteration that populates global_map_u / global_map_p
      // in CumulantSolver::solve, so G0n receives operators in the same order the old path used.
      std::vector<int> u_global_idx;
      std::vector<int> p_global_idx;
    };

    struct ProductTerm {
      int sign;                            // ±1, already = (-1) * sign_u * sign_p
      std::vector<int> factor_node_ids;    // each id < index of the parent node (topological order)
    };

    struct Node {
      LeafOps leaf;
      std::vector<ProductTerm> subtraction_terms;
    };

    std::vector<Node> nodes;   // topologically ordered: children before parent
    int root_id = -1;          // always nodes.size() - 1 after a successful build
  };

  // Evaluate a plan against the current τ values and operators of the vertex.
  // master_unprimed.taus[u_global_idx[j]] gives the τ for leaf slot j; same for ops.
  // solver + infinite_U select which G0n variant to call at each leaf.
  template <int N_sites, typename T>
  T evaluate_plan(CumulantPlan const &plan,
                  Args<N_sites, T> const &master_unprimed,
                  Args<N_sites, T> const &master_primed,
                  HubbardSolver<N_sites, T> const &solver,
                  bool infinite_U);

} // namespace sc_expansion
