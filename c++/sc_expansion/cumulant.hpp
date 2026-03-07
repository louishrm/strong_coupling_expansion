#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_atom.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace sc_expansion {

  template <typename T>
  class CumulantSolver {
    public:
    using Arg     = std::pair<double, int>;
    using ArgList = std::vector<Arg>;

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
    const HubbardAtom<T> &atom;
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
    CumulantSolver(const ArgList &u, const ArgList &p, const HubbardAtom<T> &a, bool infinite_U);

    mutable int cache_hits   = 0;
    mutable int cache_misses = 0;

    // --- The Core Recursive Function ---
    T solve(uint64_t mask_u, uint64_t mask_p);

    // Wrapper to match the test usage
    T compute_cumulant_decomposition();
  };

  template <typename T>
  T compute_cumulant_decomposition(ArgList const &unprimed, ArgList const &primed, HubbardAtom<T> const &atom,
                                        bool infinite_U = false, bool verbose = false);

  // Wrapper class for Python to easily compute cumulants for spin 0
  class CumulantHelper {
    HubbardAtom<double> atom;

    public:
    CumulantHelper(double U, double beta, double mu) : atom(U, beta, mu) {}

    double compute(std::vector<double> taus) {
      if (taus.size() % 2 != 0) { throw std::invalid_argument("CumulantHelper::compute: input vector must have even size (2n)."); }
      size_t n = taus.size() / 2;

      ArgList u, p;
      u.reserve(n);
      p.reserve(n);

      for (size_t i = 0; i < n; ++i) {
        u.push_back({taus[i], 0});     // First n elements
        p.push_back({taus[n + i], 0}); // Last n elements
      }
      return compute_cumulant_decomposition(u, p, atom);
    }
  };

} // namespace sc_expansion

#endif