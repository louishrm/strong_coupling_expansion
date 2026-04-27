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
      std::vector<ProductTerm> subtraction_terms;
    };

    std::vector<Node> nodes;
    int root_id = -1;
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

    private:
    struct CacheKey {
      uint64_t u_mask;
      uint64_t p_mask;
      bool operator==(CacheKey const &o) const { return this->u_mask == o.u_mask && this->p_mask == o.p_mask; }
    };

    struct KeyHasher {
      std::size_t operator()(CacheKey const &k) const {
        return std::hash<uint64_t>{}(k.u_mask) ^ (std::hash<uint64_t>{}(k.p_mask) << 1);
      }
    };

    Args<N_sites, T> const &master_unprimed;
    Args<N_sites, T> const &master_primed;
    HubbardSolver<N_sites, T> const *solver_ptr = nullptr;
    bool infinite_U                             = false;

    std::unordered_map<CacheKey, int, KeyHasher> plan_node_ids;
    uint64_t plan_spin_mask_stable_u = 0;
    uint64_t plan_spin_mask_stable_p = 0;
    std::vector<int> plan_inv_argsort_u;
    std::vector<int> plan_inv_argsort_p;

    int solve_record(uint64_t mask_u, uint64_t mask_p, CumulantPlan &plan);

    void record_distribute_primed(std::vector<uint64_t> const &u_partition_masks, int u_idx, uint64_t current_p_pool, int overall_sign,
                                  std::vector<int> &factors_so_far, std::vector<int> const &stable_map_u, std::vector<int> const &stable_map_p,
                                  std::vector<CumulantPlan::ProductTerm> &out, CumulantPlan &plan);
  };

  template <int N_sites, typename T>
  T evaluate_plan(CumulantPlan const &plan, Args<N_sites, T> const &master_unprimed, Args<N_sites, T> const &master_primed,
                  HubbardSolver<N_sites, T> const &solver, bool infinite_U);

} // namespace sc_expansion
