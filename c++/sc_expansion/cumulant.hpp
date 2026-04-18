#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_solver.hpp"
#include "args.hpp"
#include "cumulant_plan.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace sc_expansion {

  template <int N_sites, typename T>
  class CumulantSolver {
    public:
    // Cache Key: pair<unprimed_mask, primed_mask>
    struct CacheKey {
      uint64_t u_mask;
      uint64_t p_mask;

      bool operator==(const CacheKey &o) const { return this->u_mask == o.u_mask && this->p_mask == o.p_mask; }
    };

    struct KeyHasher {
      std::size_t operator()(const CacheKey &k) const {
        // simple hash combine
        return std::hash<uint64_t>{}(k.u_mask) ^ (std::hash<uint64_t>{}(k.p_mask) << 1);
      }
    };

    private:
    // References to the original full lists (The "Master" lists)
    const Args<N_sites, T> &master_unprimed;
    const Args<N_sites, T> &master_primed;
    const HubbardSolver<N_sites, T> &solver;
    bool infinite_U = false;

    T call_bare(const Args<N_sites, T> &u, const Args<N_sites, T> &p) const;

    // Memoization Table
    std::unordered_map<CacheKey, T, KeyHasher> memo;

    // Pre-calculated spin masks for fast conservation checks
    uint64_t master_spin_mask_u = 0;
    uint64_t master_spin_mask_p = 0;

    T distribute_primed(const std::vector<uint64_t> &u_partition_masks, int u_idx, uint64_t current_p_pool, const std::vector<int> &global_map_u,
                             const std::vector<int> &global_map_p);

    public:
    CumulantSolver(const Args<N_sites, T> &u, const Args<N_sites, T> &p, const HubbardSolver<N_sites, T> &s, bool infinite_U);

    mutable int cache_hits   = 0;
    mutable int cache_misses = 0;

    // --- The Core Recursive Function ---
    T solve(uint64_t mask_u, uint64_t mask_p);

    // Wrapper to match the test usage
    T compute_cumulant_decomposition();

    // --- Plan recording (τ-independent decomposition blueprint) ---
    // Walks the same recursion as compute_cumulant_decomposition but, instead of evaluating
    // numerical cumulant values, records the decomposition structure into `plan`.
    // After this call, evaluate_plan(plan, ...) reproduces the numerical result at any τ.
    void record_plan(CumulantPlan &plan);

    private:
    // Memo mapping (mask_u, mask_p) → node id in the plan being built.
    // The masks are over STABLE (pre-sort) positions — see note on record_plan.
    // Distinct from `memo` (which stores T values keyed on sorted-basis masks). -1 = zero.
    std::unordered_map<CacheKey, int, KeyHasher> plan_node_ids;

    // Stable-basis spin masks and stable→sorted permutation, populated at the start of
    // record_plan(). solve_record / record_distribute_primed use STABLE masks internally
    // so that all combinatorial signs are τ-invariant.
    uint64_t plan_spin_mask_stable_u = 0;
    uint64_t plan_spin_mask_stable_p = 0;
    std::vector<int> plan_inv_argsort_u; // inv_argsort_u[stable_pos] = sorted_pos in master_unprimed
    std::vector<int> plan_inv_argsort_p;

    // Mirrors solve(). Returns the node id representing C(mask_u|mask_p), or -1 if zero.
    int solve_record(uint64_t mask_u, uint64_t mask_p, CumulantPlan &plan);

    // Mirrors distribute_primed(). For one fixed unprimed partition, enumerates every
    // primed-assignment and appends one ProductTerm per assignment into `out`.
    // overall_sign is the cumulative (sign_u * Π step_sign_p) carried through the recursion;
    // factors_so_far is the list of child node ids accumulated so far.
    // Returns false if any emitted term would be zero (current path is pruned).
    void record_distribute_primed(const std::vector<uint64_t> &u_partition_masks, int u_idx,
                                  uint64_t current_p_pool, int overall_sign,
                                  std::vector<int> &factors_so_far,
                                  const std::vector<int> &global_map_u,
                                  const std::vector<int> &global_map_p,
                                  std::vector<CumulantPlan::ProductTerm> &out,
                                  CumulantPlan &plan);
  };

  template <int N_sites, typename T>
  T compute_cumulant_decomposition(Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed, HubbardSolver<N_sites, T> const &solver,
                                        bool infinite_U = false, bool verbose = false);

} // namespace sc_expansion

#endif
