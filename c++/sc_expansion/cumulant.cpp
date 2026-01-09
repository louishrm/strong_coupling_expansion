#include "cumulant.hpp"
#include <numeric>
#include <set>
#include <algorithm>
#include <stdexcept>

namespace {

  // Type definitions for clarity
  using subset_t         = std::vector<int>;
  using partition_t      = std::vector<subset_t>;
  using all_partitions_t = std::vector<partition_t>;

  void fill_partitions(std::vector<int> const &set, int index, all_partitions_t &ans, partition_t &current_partition) {

    // Base Case: All elements processed
    if (index == set.size()) {
      ans.push_back(current_partition);
      return;
    }

    // Option 1: Add the current element to an existing subset
    for (size_t i = 0; i < current_partition.size(); ++i) {
      current_partition[i].push_back(set[index]);              // Do
      fill_partitions(set, index + 1, ans, current_partition); // Recurse
      current_partition[i].pop_back();                         // Undo (Backtrack)
    }

    // Option 2: Add the current element as a new singleton subset
    current_partition.push_back({set[index]});               // Do
    fill_partitions(set, index + 1, ans, current_partition); // Recurse
    current_partition.pop_back();                            // Undo (Backtrack)
  }

  all_partitions_t get_partitions(std::vector<int> const &set) {

    all_partitions_t ans;
    partition_t current_partition;
    fill_partitions(set, 0, ans, current_partition);
    return ans;
  }

  std::vector<std::vector<int>> generate_permutations(std::vector<int> &set) {

    std::vector<std::vector<int>> permutations;
    std::sort(set.begin(), set.end());
    do { permutations.push_back(set); } while (std::next_permutation(set.begin(), set.end()));
    return permutations;
  }

  int permutation_sign(std::vector<int> const &perm) {

    int inversions = 0;
    for (size_t i = 0; i < perm.size(); ++i) {
      for (size_t j = i + 1; j < perm.size(); ++j) {
        if (perm[i] > perm[j]) { ++inversions; }
      }
    }
    return (inversions % 2 == 0) ? 1 : -1;
  }

  bool is_spin_conserving(sc_expansion::HubbardAtom::cumul_args const &unprimed, sc_expansion::HubbardAtom::cumul_args const &primed) {
    int up_count_unprimed   = 0;
    int down_count_unprimed = 0;
    for (const auto &arg : unprimed) {
      if (arg.second == 0)
        up_count_unprimed++;
      else
        down_count_unprimed++;
    }

    int up_count_primed   = 0;
    int down_count_primed = 0;
    for (const auto &arg : primed) {
      if (arg.second == 0)
        up_count_primed++;
      else
        down_count_primed++;
    }

    return (up_count_unprimed == up_count_primed) && (down_count_unprimed == down_count_primed);
  }

} // namespace

namespace sc_expansion {

  CumulantSolver::CumulantSolver(const ArgList &u, const ArgList &p, const HubbardAtom &a) : master_unprimed(u), master_primed(p), atom(a) {
    // Pre-calculate spin patterns for fast checks
    // Assuming Arg.second is the spin (0 or 1)
    for (size_t i = 0; i < u.size(); ++i) {
      if (u[i].second == 1) master_spin_mask_u |= (1ULL << i);
    }
    for (size_t i = 0; i < p.size(); ++i) {
      if (p[i].second == 1) master_spin_mask_p |= (1ULL << i);
    }
  }

  double CumulantSolver::solve(uint64_t mask_u, uint64_t mask_p) {

    // 1. Check Spin Conservation (Fast Bitwise Check)
    int up_u = __builtin_popcountll(mask_u & master_spin_mask_u);
    int up_p = __builtin_popcountll(mask_p & master_spin_mask_p);
    if (up_u != up_p) return 0.0;

    // 2. Check Cache
    // 2. Check Cache
    CacheKey key{mask_u, mask_p};
    if (auto it = memo.find(key); it != memo.end()) {
      cache_hits++; // <--- HIT!
      return it->second;
    }

    cache_misses++; // <--- MISS (We have to compute it)

    // 3. Prepare Mapping: Global Indices -> Local Indices (0..k-1)
    std::vector<int> global_map_u;
    std::vector<int> global_map_p;

    // Extract set bits to create the mapping
    for (int i = 0; i < 64; ++i) {
      if ((mask_u >> i) & 1) global_map_u.push_back(i);
    }
    for (int i = 0; i < 64; ++i) {
      if ((mask_p >> i) & 1) global_map_p.push_back(i);
    }

    // 4. Base Case: Order 1
    if (global_map_u.size() == 1) {
      // Reconstruct temp vectors for the atom call
      ArgList args_u   = {master_unprimed[global_map_u[0]]};
      ArgList args_p   = {master_primed[global_map_p[0]]};
      return memo[key] = atom.G0(args_u, args_p);
    }

    // 5. Compute G0_n (First term)
    ArgList current_args_u, current_args_p;
    current_args_u.reserve(global_map_u.size());
    current_args_p.reserve(global_map_p.size());

    for (int idx : global_map_u) current_args_u.push_back(master_unprimed[idx]);
    for (int idx : global_map_p) current_args_p.push_back(master_primed[idx]);

    double G0n = atom.G0(current_args_u, current_args_p);

    // 6. Subtraction Term Logic
    double low_order_cumulants = 0.0;
    int order                  = global_map_u.size();

    // Create local indices 0, 1, ..., k-1
    std::vector<int> local_indices(order);
    std::iota(local_indices.begin(), local_indices.end(), 0);

    auto unprimed_partitions = get_partitions(local_indices);
    auto primed_perms        = generate_permutations(local_indices);

    for (const auto &current_partition : unprimed_partitions) {
      if (current_partition.size() == 1) continue;

      std::vector<int> flattened_unprimed_perm;
      flattened_unprimed_perm.reserve(order);

      for (auto const &subset : current_partition) { flattened_unprimed_perm.insert(flattened_unprimed_perm.end(), subset.begin(), subset.end()); }

      int sign_u = permutation_sign(flattened_unprimed_perm);

      std::set<std::vector<int>> visited_permutations;

      for (const auto &perm : primed_perms) {
        int primed_cursor = 0;
        std::vector<int> effective_permutation(order);
        std::vector<std::vector<int>> primed_subsets;
        auto write_cursor = effective_permutation.begin();

        for (const auto &x : current_partition) {
          auto write_end_cursor = write_cursor + x.size();
          std::copy(perm.begin() + primed_cursor, perm.begin() + primed_cursor + x.size(), write_cursor);
          std::sort(write_cursor, write_end_cursor);
          primed_subsets.push_back(std::vector<int>(write_cursor, write_end_cursor));
          primed_cursor += x.size();
          write_cursor = write_end_cursor;
        }

        if (visited_permutations.count(effective_permutation)) continue;
        visited_permutations.insert(effective_permutation);

        double current_product = 1.0;
        int sign_p             = permutation_sign(effective_permutation);
        int sign               = -sign_p * sign_u;

        for (size_t i = 0; i < current_partition.size(); ++i) {
          const auto &local_u_subset = current_partition[i];
          const auto &local_p_subset = primed_subsets[i];

          uint64_t next_mask_u = 0;
          uint64_t next_mask_p = 0;

          for (int local_idx : local_u_subset) {
            int global_idx = global_map_u[local_idx];
            next_mask_u |= (1ULL << global_idx);
          }

          for (int local_idx : local_p_subset) {
            int global_idx = global_map_p[local_idx];
            next_mask_p |= (1ULL << global_idx);
          }

          current_product *= solve(next_mask_u, next_mask_p);
        }
        low_order_cumulants += sign * current_product;
      }
    }

    return memo[key] = G0n + low_order_cumulants;
  }

  double CumulantSolver::compute_cumulant_decomposition() {
    uint64_t full_mask = (master_unprimed.size() == 64) ? ~0ULL : (1ULL << master_unprimed.size()) - 1;
    return solve(full_mask, full_mask);
  }

  double compute_cumulant_decomposition(HubbardAtom::cumul_args const &unprimed, HubbardAtom::cumul_args const &primed, HubbardAtom const &atom,
                                        bool verbose) {
    if (unprimed.size() != primed.size())
      throw std::invalid_argument("CumulantSolver::compute_cumulant_decomposition: unprimed and primed lists must have the same size.");
    if (unprimed.empty()) throw std::invalid_argument("CumulantSolver::compute_cumulant_decomposition: unprimed and primed lists cannot be empty.");

    CumulantSolver solver(unprimed, primed, atom);

    double result = solver.compute_cumulant_decomposition();

    if (verbose) { std::cout << "CumulantSolver Cache Stats: Hits = " << solver.cache_hits << ", Misses = " << solver.cache_misses << "\n"; }

    return result;
  }

} // namespace sc_expansion
