#include "cumulant.hpp"
#include <numeric>
#include <set>
#include <algorithm>
#include <stdexcept>
#include <functional> // Required for std::forward

namespace {

  // Type definitions for partition generation
  using subset_t         = std::vector<int>;
  using partition_t      = std::vector<subset_t>;
  using all_partitions_t = std::vector<partition_t>;

  void fill_partitions(std::vector<int> const &set, int index, all_partitions_t &ans, partition_t &current_partition) {
    if (index == set.size()) {
      ans.push_back(current_partition);
      return;
    }
    for (size_t i = 0; i < current_partition.size(); ++i) {
      current_partition[i].push_back(set[index]);
      fill_partitions(set, index + 1, ans, current_partition);
      current_partition[i].pop_back();
    }
    current_partition.push_back({set[index]});
    fill_partitions(set, index + 1, ans, current_partition);
    current_partition.pop_back();
  }

  const all_partitions_t &get_partitions(int n) {
    static std::unordered_map<int, all_partitions_t> partition_cache;
    auto it = partition_cache.find(n);
    if (it != partition_cache.end()) return it->second;

    std::vector<int> set(n);
    std::iota(set.begin(), set.end(), 0);
    all_partitions_t ans;
    partition_t current_partition;
    fill_partitions(set, 0, ans, current_partition);
    return partition_cache[n] = std::move(ans);
  }

  // --- Bitwise Helpers ---

  inline int popcount(uint64_t n) { return __builtin_popcountll(n); }
  inline int ctz(uint64_t n) { return __builtin_ctzll(n); }

  int compute_extraction_sign(uint64_t pool, uint64_t subset) {
    int inversions     = 0;
    uint64_t remaining = pool ^ subset;

    while (subset) {
      int bit             = ctz(subset);
      uint64_t lower_mask = (1ULL << bit) - 1;
      inversions += popcount(remaining & lower_mask);
      subset &= ~(1ULL << bit);
    }
    return (inversions % 2 == 0) ? 1 : -1;
  }

  // Recursive Subset Generator (Avoids template depth issues)
  template <typename Func> void recursive_subset_generator(uint64_t pool, int k, uint64_t current_subset, Func &&callback) {
    if (k == 0) {
      callback(current_subset);
      return;
    }
    if (popcount(pool) < k) return;

    int bit            = ctz(pool);
    uint64_t bit_mask  = (1ULL << bit);
    uint64_t next_pool = pool ^ bit_mask;

    // Branch A: Take bit
    recursive_subset_generator(next_pool, k - 1, current_subset | bit_mask, callback);

    // Branch B: Skip bit
    recursive_subset_generator(next_pool, k, current_subset, callback);
  }

  // Public wrapper
  template <typename Func> void for_each_subset(uint64_t pool, int k, Func &&callback) {
    recursive_subset_generator(pool, k, 0, std::forward<Func>(callback));
  }

} // namespace

namespace sc_expansion {

  CumulantSolver::CumulantSolver(const ArgList &u, const ArgList &p, const HubbardAtom &a, bool infinite_U_)
     : master_unprimed(u), master_primed(p), atom(a), infinite_U(infinite_U_) {

    // Pre-calculate spin masks
    for (size_t i = 0; i < u.size(); ++i) {
      if (u[i].second == 1) this->master_spin_mask_u |= (1ULL << i);
    }
    for (size_t i = 0; i < p.size(); ++i) {
      if (p[i].second == 1) this->master_spin_mask_p |= (1ULL << i);
    }
  }

  double CumulantSolver::call_bare(const ArgList &u, const ArgList &p) const {
    int n = u.size();
    std::vector<double> taus(2 * n);
    std::vector<int> spins(2 * n);
    for (int i = 0; i < n; ++i) {
      // Creation (primed) goes to even indices
      taus[2 * i]  = p[i].first;
      spins[2 * i] = p[i].second;
      // Destruction (unprimed) goes to odd indices
      taus[2 * i + 1]  = u[i].first;
      spins[2 * i + 1] = u[i].second;
    }
    return infinite_U ? atom.G0_infinite_U(taus, spins) : atom.G0(taus, spins);
  }

  // Recursive Distributor (Fixed Sign Logic)
  double CumulantSolver::distribute_primed(const std::vector<uint64_t> &u_partition_masks, int u_idx, uint64_t current_p_pool,
                                           const std::vector<int> &global_map_u, const std::vector<int> &global_map_p) {

    // Base Case
    if (u_idx == u_partition_masks.size()) { return 1.0; }

    uint64_t u_mask  = u_partition_masks[u_idx];
    int needed_k     = popcount(u_mask);
    double sum_terms = 0.0;

    for_each_subset(current_p_pool, needed_k, [&](uint64_t p_submask) {
      // 1. Calculate Local Sign for this step
      int step_sign_p = compute_extraction_sign(current_p_pool, p_submask);

      // 2. Map Local -> Global
      uint64_t global_mask_u = 0;
      uint64_t global_mask_p = 0;

      uint64_t temp_u = u_mask;
      while (temp_u) {
        int b = ctz(temp_u);
        global_mask_u |= (1ULL << global_map_u[b]);
        temp_u &= ~(1ULL << b);
      }

      uint64_t temp_p = p_submask;
      while (temp_p) {
        int b = ctz(temp_p);
        global_mask_p |= (1ULL << global_map_p[b]);
        temp_p &= ~(1ULL << b);
      }

      // 3. Compute Value
      double term_val = this->solve(global_mask_u, global_mask_p);

      // 4. Recurse (No sign passed down)
      double remainder = distribute_primed(u_partition_masks, u_idx + 1, current_p_pool ^ p_submask, global_map_u, global_map_p);

      // 5. Accumulate: Sign * Value * Rest
      sum_terms += step_sign_p * term_val * remainder;
    });

    return sum_terms;
  }

  double CumulantSolver::solve(uint64_t mask_u, uint64_t mask_p) {

    // 1. Check Spin Conservation
    if (popcount(mask_u & master_spin_mask_u) != popcount(mask_p & master_spin_mask_p)) return 0.0;

    // 2. Cache Check
    CacheKey key{mask_u, mask_p};
    if (auto it = memo.find(key); it != memo.end()) {
      cache_hits++;
      return it->second;
    }
    cache_misses++;

    // 3. Map Global -> Local
    std::vector<int> global_map_u, global_map_p;
    uint64_t temp = mask_u;
    while (temp) {
      int i = ctz(temp);
      global_map_u.push_back(i);
      temp &= ~(1ULL << i);
    }
    temp = mask_p;
    while (temp) {
      int i = ctz(temp);
      global_map_p.push_back(i);
      temp &= ~(1ULL << i);
    }

    int order = global_map_u.size();

    // 4. Base Case: Order 1
    if (order == 1) {
      ArgList args_u   = {master_unprimed[global_map_u[0]]};
      ArgList args_p   = {master_primed[global_map_p[0]]};
      return memo[key] = this->call_bare(args_u, args_p);
    }

    // 5. Compute G0_n
    ArgList current_args_u, current_args_p;
    current_args_u.reserve(order);
    current_args_p.reserve(order);
    for (int idx : global_map_u) current_args_u.push_back(master_unprimed[idx]);
    for (int idx : global_map_p) current_args_p.push_back(master_primed[idx]);

    double G0n = this->call_bare(current_args_u, current_args_p);

    // 6. Subtraction Term Logic
    double low_order_cumulants = 0.0;
    const auto &unprimed_partitions = get_partitions(order);

    for (const auto &partition : unprimed_partitions) {
      if (partition.size() == 1) continue;

      int sign_u      = 1;
      uint64_t u_pool = (order == 64) ? ~0ULL : (1ULL << order) - 1;
      std::vector<uint64_t> u_partition_masks;

      // Calculate S_u (Unprimed Sign)
      for (const auto &subset : partition) {
        uint64_t submask = 0;
        for (int idx : subset) submask |= (1ULL << idx);

        sign_u *= compute_extraction_sign(u_pool, submask);
        u_pool ^= submask;

        u_partition_masks.push_back(submask);
      }

      // Calculate Sum of (S_p * Terms)
      uint64_t full_local_p_mask = (order == 64) ? ~0ULL : (1ULL << order) - 1;

      double term_sum = distribute_primed(u_partition_masks, 0, full_local_p_mask, global_map_u, global_map_p);

      // Apply Global Sign
      low_order_cumulants += (-1.0 * sign_u) * term_sum;
    }

    return memo[key] = G0n + low_order_cumulants;
  }

  double CumulantSolver::compute_cumulant_decomposition() {
    uint64_t full_mask = (this->master_unprimed.size() == 64) ? ~0ULL : (1ULL << this->master_unprimed.size()) - 1;
    return this->solve(full_mask, full_mask);
  }

  double compute_cumulant_decomposition(ArgList const &unprimed, ArgList const &primed, HubbardAtom const &atom,
                                        bool infinite_U, bool verbose) {
    if (unprimed.size() != primed.size()) throw std::invalid_argument("Size mismatch in compute_cumulant_decomposition");
    if (unprimed.empty()) throw std::invalid_argument("Empty list in compute_cumulant_decomposition");

    CumulantSolver solver(unprimed, primed, atom, infinite_U);
    double result = solver.compute_cumulant_decomposition();

    if (verbose) { std::cout << "CumulantSolver Cache Stats: Hits = " << solver.cache_hits << ", Misses = " << solver.cache_misses << "\n"; }
    return result;
  }

} // namespace sc_expansion