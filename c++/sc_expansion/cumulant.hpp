#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_solver.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace sc_expansion {

  using Arg     = std::pair<double, int>;
  using ArgList = std::vector<Arg>;

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
    const ArgList &master_unprimed;
    const ArgList &master_primed;
    const HubbardSolver<N_sites, T> &solver;
    bool infinite_U = false;

    T call_bare(const ArgList &u, const ArgList &p) const;

    // Memoization Table
    std::unordered_map<CacheKey, T, KeyHasher> memo;

    // Pre-calculated spin masks for fast conservation checks
    uint64_t master_spin_mask_u = 0;
    uint64_t master_spin_mask_p = 0;

    T distribute_primed(const std::vector<uint64_t> &u_partition_masks, int u_idx, uint64_t current_p_pool, const std::vector<int> &global_map_u,
                             const std::vector<int> &global_map_p);

    public:
    CumulantSolver(const ArgList &u, const ArgList &p, const HubbardSolver<N_sites, T> &s, bool infinite_U);

    mutable int cache_hits   = 0;
    mutable int cache_misses = 0;

    // --- The Core Recursive Function ---
    T solve(uint64_t mask_u, uint64_t mask_p);

    // Wrapper to match the test usage
    T compute_cumulant_decomposition();
  };

  template <int N_sites, typename T>
  T compute_cumulant_decomposition(ArgList const &unprimed, ArgList const &primed, HubbardSolver<N_sites, T> const &solver,
                                        bool infinite_U = false, bool verbose = false);

} // namespace sc_expansion

#endif
