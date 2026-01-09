#ifndef CUMULANT_HPP
#define CUMULANT_HPP

#include "hubbard_atom.hpp"
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace sc_expansion {

  class CumulantSolver {
    public:
    using Arg     = std::pair<double, int>;
    using ArgList = std::vector<Arg>;

    // Cache Key: pair<unprimed_mask, primed_mask>
    struct CacheKey {
      uint64_t u_mask;
      uint64_t p_mask;

      bool operator==(const CacheKey &o) const { return u_mask == o.u_mask && p_mask == o.p_mask; }
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
    const HubbardAtom &atom;

    // Memoization Table
    std::unordered_map<CacheKey, double, KeyHasher> memo;

    // Pre-calculated spin masks for fast conservation checks
    uint64_t master_spin_mask_u = 0;
    uint64_t master_spin_mask_p = 0;

    public:
    CumulantSolver(const ArgList &u, const ArgList &p, const HubbardAtom &a);

    mutable int cache_hits   = 0;
    mutable int cache_misses = 0;

    // --- The Core Recursive Function ---
    double solve(uint64_t mask_u, uint64_t mask_p);

    // Wrapper to match the test usage
    double compute_cumulant_decomposition();
  };

  double compute_cumulant_decomposition(HubbardAtom::cumul_args const &unprimed, HubbardAtom::cumul_args const &primed, HubbardAtom const &atom,
                                        bool verbose);

} // namespace sc_expansion

#endif