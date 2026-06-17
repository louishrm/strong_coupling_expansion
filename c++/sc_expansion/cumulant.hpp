#pragma once

#include "hubbard_solver.hpp"
#include "args.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace sc_expansion {

  // τ-independent decomposition blueprint for one vertex.
  // `nodes[root_id]` yields the full cumulant: G0n at the root leaf minus a sum of products
  // of lower-order cumulants (subtraction_terms), each a signed product of child node values.
  struct CumulantPlan {

    struct LeafOps {
      std::vector<int> u_global_idx;
      std::vector<int> p_global_idx;
    };

    struct ProductTerm {
      int sign;
      std::vector<int> factor_node_ids;
    };

    struct Node {
      LeafOps leaf;
      // Static-density decoration attached to this sub-cumulant's leaf.
      // Empty  → ordinary G0n leaf (existing behavior).
      // Non-empty → G0n_with_densities leaf: each entry is one bosonic
      // density atom n_σ(0) bound to this subset. Multi-entry is allowed
      // (e.g. two same-vertex marks co-routed into one block). The atom is
      // NOT a hybridization leg; it does not enter the partition lattice.
      std::vector<int> leaf_density_orbitals;
      std::vector<ProductTerm> subtraction_terms;
    };

    std::vector<Node> nodes;
    int root_id = -1;

    // For rooted (partial-cumulant) plans built with one or more block
    // constraints: each size-1 sub-cumulant whose only ops are one mark's
    // (c†_σ, c_σ) pair has its leaf value replaced by ⟨n_σ⟩ at evaluation
    // time (bypasses the equal-time T_τ ambiguity that ⟨T_τ c_σ(0) c†_σ(0)⟩
    // would otherwise introduce). One entry per block.
    std::vector<int> mark_node_ids;
    std::vector<int> mark_orbitals;

    // Coincidence groups: sets of mark blocks (registered via
    // add_coincidence_group) whose insertions sit at the same physical
    // location and spin. For same-spin coincident marks the product
    // n_σ(0)^k = n_σ(0) collapses; in mixed leaves we drop k-1 redundant
    // (c†_σ, c_σ) pairs before invoking G0n. Pure-group leaves of size
    // >= 2 are handled via the mark_node_ids/mark_orbitals override.
    struct CoincidenceGroup {
      int orbital;
      std::vector<int> u_stable; // stable indices of the group's unprimed legs
      std::vector<int> p_stable; // stable indices of the group's primed legs
    };
    std::vector<CoincidenceGroup> coincidence_groups;
  };

  // Records a CumulantPlan for a given (unprimed, primed) vertex leg structure.
  // The plan is independent of τ values and of infinite_U, so it can be reused across
  // all Monte Carlo steps and both finite/infinite-U evaluations.
  template <int N_sites, typename T> class CumulantSolver {
    public:
    CumulantSolver(Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed);

    // Solver-aware overload: also retains a HubbardSolver pointer + infinite_U
    // flag so callers can ask for the cumulant value directly via
    // compute_cumulant_decomposition() without manually managing a CumulantPlan.
    CumulantSolver(Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed, HubbardSolver<N_sites, T> const &solver, bool infinite_U);

    void record_plan(CumulantPlan &plan);

    // Build a plan and evaluate it in one shot, using the solver+infinite_U
    // captured at construction. Intended as a reference path for tests; hot
    // MC code should call record_plan once and evaluate_plan() per step.
    T compute_cumulant_decomposition();

    // Pin a pair of operators (one unprimed, one primed) into the same
    // partition block. Used for rooted (correlator) vertices, where each
    // mark insertion n_σ(0) = c†_σ(0)·c_σ(0) acts as one composite operator
    // in the cumulant partition lattice — producing a partial cumulant in
    // place of the full κ_n. Call once per mark; multiple blocks can coexist
    // on the same vertex (same-site density-density insertions).
    // Indices are *input* positions (pre-sort) within the unprimed / primed
    // operator lists.
    void add_block_constraint(int u_input_idx, int p_input_idx) {
      this->block_u_stable_list.push_back(u_input_idx);
      this->block_p_stable_list.push_back(p_input_idx);
    }

    // Declare a group of coincident mark blocks (block indices reference
    // the order of prior add_block_constraint calls). All marks in a group
    // sit at the same (τ=0, σ=orbital) — their product n_σ^k collapses to
    // n_σ. Recorded into the plan for use at evaluation.
    void add_coincidence_group(int orbital, std::vector<int> block_indices) {
      this->pending_coincidence_groups.push_back({orbital, std::move(block_indices)});
    }

    // Register one static density n_σ(0) decoration on this vertex. Unlike
    // add_block_constraint (which folds the density's c†/c into the
    // operator lists as a paired hybridization block), this records the
    // density as an external decoration: the underlying c†_σ(0), c_σ(0)
    // never enter master_unprimed / master_primed; instead the orbital is
    // attached to the cumulant plan's leaves via leaf_density_orbitals.
    // Multiple calls accumulate (e.g. for two coincident marks on the
    // same vertex). Order does not matter.
    void add_static_density(int orbital) { this->pending_static_densities.push_back(orbital); }

    private:
    struct CacheKey {
      uint64_t u_mask;
      uint64_t p_mask;
      // Compact multiset encoding of the sub-cumulant's density decoration:
      // 4 bits per orbital index (count of n_σ insertions on this sub-cumulant).
      // 0 = no decoration (matches the original behavior).
      uint64_t decoration_key;
      bool operator==(CacheKey const &o) const {
        return this->u_mask == o.u_mask && this->p_mask == o.p_mask && this->decoration_key == o.decoration_key;
      }
    };

    struct KeyHasher {
      std::size_t operator()(CacheKey const &k) const {
        return std::hash<uint64_t>{}(k.u_mask) ^ (std::hash<uint64_t>{}(k.p_mask) << 1) ^ (std::hash<uint64_t>{}(k.decoration_key) << 2);
      }
    };

    Args<N_sites, T> const &master_unprimed;
    Args<N_sites, T> const &master_primed;
    HubbardSolver<N_sites, T> const *solver_ptr = nullptr;
    bool infinite_U                             = false;

    std::unordered_map<CacheKey, int, KeyHasher> plan_node_ids;
    uint64_t plan_spin_mask_stable_u = 0;
    uint64_t plan_spin_mask_stable_p = 0;

    // Block-constraint state for rooted/partial-cumulant evaluation.
    // Empty = inactive. Each entry pairs an unprimed and primed input
    // position that MUST sit together in the same partition subset; ops
    // from different blocks are free to share or split subsets independently.
    std::vector<int> block_u_stable_list;
    std::vector<int> block_p_stable_list;

    struct PendingCoincidenceGroup {
      int orbital;
      std::vector<int> block_indices;
    };
    std::vector<PendingCoincidenceGroup> pending_coincidence_groups;

    // Static density decorations registered via add_static_density(). Each
    // entry is one orbital index n_σ(0). Consumed by record_plan (Step 4)
    // to drive per-atom routing through the partition lattice.
    std::vector<int> pending_static_densities;

    std::vector<int> plan_inv_argsort_u;
    std::vector<int> plan_inv_argsort_p;

    int solve_record(uint64_t mask_u, uint64_t mask_p, std::vector<int> const &decoration, CumulantPlan &plan);

    // Collects (overall_sign, op-block-list) tuples for a fixed u-partition.
    // The op-block list is the sequence of (u_global_mask, p_global_mask) pairs
    // for each block of the joint operator partition. solve_record then loops
    // over these and enumerates per-density-atom routings on top.
    using OpBlock     = std::pair<uint64_t, uint64_t>;
    using OpBlockList = std::vector<OpBlock>;
    void record_distribute_primed_collect(std::vector<uint64_t> const &u_partition_masks, int u_idx, uint64_t current_p_pool, int overall_sign,
                                          OpBlockList &block_list_so_far, std::vector<int> const &stable_map_u, std::vector<int> const &stable_map_p,
                                          std::vector<std::pair<int, OpBlockList>> &out, std::vector<int> const &special_u_idxs,
                                          std::vector<int> const &local_p_block_idxs);

    // Recursive enumeration of how the `decoration` density atoms route across
    // the K op-blocks (existing) plus an open-ended sequence of "free groups"
    // (new leg-less blocks made of atoms only). For each routing, builds the
    // ProductTerm (calling solve_record per block) and pushes onto out_terms
    // unless the resulting joint partition has only one block (= the leaf,
    // not a subtraction).
    void enumerate_atom_routings_and_emit(std::vector<int> const &decoration, int atom_idx, OpBlockList const &op_blocks, int base_sign,
                                          std::vector<std::vector<int>> &per_op_block_decoration, std::vector<std::vector<int>> &free_groups,
                                          std::vector<CumulantPlan::ProductTerm> &out_terms, CumulantPlan &plan);
  };

  template <int N_sites, typename T>
  T evaluate_plan(CumulantPlan const &plan, Args<N_sites, T> const &master_unprimed, Args<N_sites, T> const &master_primed,
                  HubbardSolver<N_sites, T> const &solver, bool infinite_U);

} // namespace sc_expansion
