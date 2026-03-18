#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_solver.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace sc_expansion {

  template <typename T>
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
    // Master lists derived from taus and ops
    std::vector<double> master_taus;
    std::vector<int> master_ops;
    const HubbardAtom<T> &atom;
    bool infinite_U = false;

    T call_bare(uint64_t mask_u, uint64_t mask_p) const;

    // Memoization Table
    std::unordered_map<CacheKey, T, KeyHasher> memo;

    // Pre-calculated spin masks for fast conservation checks
    uint64_t master_spin_mask_u = 0;
    uint64_t master_spin_mask_p = 0;

    T distribute_primed(const std::vector<uint64_t> &u_partition_masks, int u_idx, uint64_t current_p_pool, const std::vector<int> &global_map_u,
                             const std::vector<int> &global_map_p);

    public:
    CumulantSolver(const std::vector<double> &taus, const std::vector<int> &ops, const HubbardAtom<T> &a, bool infinite_U);

    mutable int cache_hits   = 0;
    mutable int cache_misses = 0;

    // --- The Core Recursive Function ---
    T solve(uint64_t mask_u, uint64_t mask_p);

    // Wrapper to match the test usage
    T compute_cumulant_decomposition();
  };

  template <typename T>
  T compute_cumulant_decomposition(std::vector<double> const &taus, std::vector<int> const &ops, HubbardAtom<T> const &atom,
                                        bool infinite_U = false, bool verbose = false);

  // Wrapper class for Python to easily compute cumulants for spin 0
  class CumulantHelper {
    Parameters<double> params;
    HubbardAtom<double> atom;

    public:
    CumulantHelper(double U, double beta, double mu) : params{U, beta, mu}, atom(params) {}

    double compute(std::vector<double> taus) {
      if (taus.size() % 2 != 0) { throw std::invalid_argument("CumulantHelper::compute: input vector must have even size (2n)."); }
      size_t n = taus.size() / 2;

      std::vector<double> inter_taus(2 * n);
      std::vector<int> inter_ops(2 * n);

      for (size_t i = 0; i < n; ++i) {
        inter_taus[2 * i + 1] = taus[i]; // unprimed
        inter_ops[2 * i + 1] = 0;        // spin 0, action 0 (destroy)
        inter_taus[2 * i] = taus[n + i]; // primed
        inter_ops[2 * i] = 1;            // spin 0, action 1 (create)
      }
      return compute_cumulant_decomposition(inter_taus, inter_ops, atom);
    }
  };

} // namespace sc_expansion

#endif