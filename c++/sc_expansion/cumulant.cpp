#include "cumulant.hpp"
#include "dual.hpp"
#include <numeric>
#include <algorithm>
#include <functional>

namespace {

  using subset_t         = std::vector<int>;
  using partition_t      = std::vector<subset_t>;
  using all_partitions_t = std::vector<partition_t>;

  void fill_partitions(std::vector<int> const &set, int index, all_partitions_t &ans, partition_t &current_partition) {
    if (index == (int)set.size()) {
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

  // Enumerates all k-element subsets of `pool` (a bitmask) and invokes `callback` on each.
  template <typename Func> void recursive_subset_generator(uint64_t pool, int k, uint64_t current_subset, Func &&callback) {
    if (k == 0) {
      callback(current_subset);
      return;
    }
    if (popcount(pool) < k) return;

    int bit            = ctz(pool);
    uint64_t bit_mask  = (1ULL << bit);
    uint64_t next_pool = pool ^ bit_mask;

    recursive_subset_generator(next_pool, k - 1, current_subset | bit_mask, callback);
    recursive_subset_generator(next_pool, k, current_subset, callback);
  }

  template <typename Func> void for_each_subset(uint64_t pool, int k, Func &&callback) {
    recursive_subset_generator(pool, k, 0, std::forward<Func>(callback));
  }

} // namespace

namespace sc_expansion {

  // ============================================================================
  //  Plan recording
  // ============================================================================

  template <typename T>
  CumulantSolver<T>::CumulantSolver(Args<T> const &unprimed, Args<T> const &primed, std::vector<std::pair<int, int>> self_loop_pairs_)
     : master_unprimed(unprimed), master_primed(primed), self_loop_pairs(std::move(self_loop_pairs_)) {}

  template <typename T> uint64_t CumulantSolver<T>::forced_p_bits_for(uint64_t u_mask) const {
    uint64_t forced = 0;
    for (auto const &pair : this->self_loop_pairs) {
      if (u_mask & (1ULL << pair.first)) forced |= (1ULL << pair.second);
    }
    return forced;
  }

  template <typename T>
  void CumulantSolver<T>::record_distribute_primed(std::vector<uint64_t> const &u_partition_masks, int u_idx, uint64_t current_p_pool,
                                                   int overall_sign, std::vector<int> &factors_so_far, std::vector<int> const &stable_map_u,
                                                   std::vector<int> const &stable_map_p, std::vector<CumulantPlan::ProductTerm> &out,
                                                   CumulantPlan &plan) {

    if (u_idx == (int)u_partition_masks.size()) {
      CumulantPlan::ProductTerm term;
      term.sign            = overall_sign;
      term.factor_node_ids = factors_so_far;
      out.push_back(std::move(term));
      return;
    }

    uint64_t u_mask = u_partition_masks[u_idx];
    int needed_k    = popcount(u_mask);

    // Self-loop pairs must never be split across partition blocks: whenever a
    // u-bit belonging to a self-loop sits in `u_mask`, the matched p-bit is
    // forced into the p-submask. This gives the density-density semantics at
    // vertices that carry self-loops; at vertices without self-loops
    // `forced_p_bits` is zero and we fall back to the unrestricted enumeration.
    uint64_t forced_p_bits = this->forced_p_bits_for(u_mask);
    if ((forced_p_bits & current_p_pool) != forced_p_bits) return;
    int forced_k = popcount(forced_p_bits);
    if (forced_k > needed_k) return;
    uint64_t free_pool = current_p_pool & ~forced_p_bits;
    int free_k         = needed_k - forced_k;

    for_each_subset(free_pool, free_k, [&](uint64_t free_submask) {
      uint64_t p_submask = forced_p_bits | free_submask;
      int step_sign_p = compute_extraction_sign(current_p_pool, p_submask);

      uint64_t global_mask_u_stable = 0, global_mask_p_stable = 0;
      {
        uint64_t t = u_mask;
        while (t) { int b = ctz(t); global_mask_u_stable |= (1ULL << stable_map_u[b]); t &= ~(1ULL << b); }
      }
      {
        uint64_t t = p_submask;
        while (t) { int b = ctz(t); global_mask_p_stable |= (1ULL << stable_map_p[b]); t &= ~(1ULL << b); }
      }

      int child_id = this->solve_record(global_mask_u_stable, global_mask_p_stable, plan);
      if (child_id < 0) return;

      factors_so_far.push_back(child_id);
      this->record_distribute_primed(u_partition_masks, u_idx + 1, current_p_pool ^ p_submask, overall_sign * step_sign_p, factors_so_far,
                                     stable_map_u, stable_map_p, out, plan);
      factors_so_far.pop_back();
    });
  }

  template <typename T> int CumulantSolver<T>::solve_record(uint64_t mask_u_stable, uint64_t mask_p_stable, CumulantPlan &plan) {

    if (popcount(mask_u_stable & this->plan_spin_mask_stable_u) != popcount(mask_p_stable & this->plan_spin_mask_stable_p)) return -1;

    CacheKey key{mask_u_stable, mask_p_stable};
    if (auto it = this->plan_node_ids.find(key); it != this->plan_node_ids.end()) return it->second;

    std::vector<int> stable_map_u, stable_map_p;
    {
      uint64_t t = mask_u_stable;
      while (t) { int i = ctz(t); stable_map_u.push_back(i); t &= ~(1ULL << i); }
    }
    {
      uint64_t t = mask_p_stable;
      while (t) { int i = ctz(t); stable_map_p.push_back(i); t &= ~(1ULL << i); }
    }
    int order = (int)stable_map_u.size();

    CumulantPlan::Node node;
    node.leaf.u_global_idx = stable_map_u;
    node.leaf.p_global_idx = stable_map_p;

    if (order == 1) {
      int new_id = (int)plan.nodes.size();
      plan.nodes.push_back(std::move(node));
      this->plan_node_ids.emplace(key, new_id);
      return new_id;
    }

    auto const &unprimed_partitions = get_partitions(order);
    for (auto const &partition : unprimed_partitions) {
      if (partition.size() == 1) continue;

      int sign_u            = 1;
      uint64_t u_pool_local = (order == 64) ? ~0ULL : (1ULL << order) - 1;
      std::vector<uint64_t> u_partition_local_masks;
      for (auto const &subset : partition) {
        uint64_t submask_local = 0;
        for (int idx : subset) submask_local |= (1ULL << idx);
        sign_u *= compute_extraction_sign(u_pool_local, submask_local);
        u_pool_local ^= submask_local;
        u_partition_local_masks.push_back(submask_local);
      }

      int base_sign              = -sign_u;
      uint64_t full_local_p_mask = (order == 64) ? ~0ULL : (1ULL << order) - 1;

      std::vector<int> factors_scratch;
      factors_scratch.reserve(partition.size());
      this->record_distribute_primed(u_partition_local_masks, 0, full_local_p_mask, base_sign, factors_scratch, stable_map_u, stable_map_p,
                                     node.subtraction_terms, plan);
    }

    int new_id = (int)plan.nodes.size();
    plan.nodes.push_back(std::move(node));
    this->plan_node_ids.emplace(key, new_id);
    return new_id;
  }

  template <typename T> void CumulantSolver<T>::record_plan(CumulantPlan &plan) {
    plan.nodes.clear();
    plan.root_id = -1;
    this->plan_node_ids.clear();

    int order = this->master_unprimed.order;

    this->plan_inv_argsort_u.assign(order, 0);
    this->plan_inv_argsort_p.assign(order, 0);
    for (int sorted_pos = 0; sorted_pos < order; ++sorted_pos) {
      this->plan_inv_argsort_u[this->master_unprimed.argsort[sorted_pos]] = sorted_pos;
      this->plan_inv_argsort_p[this->master_primed.argsort[sorted_pos]]   = sorted_pos;
    }

    this->plan_spin_mask_stable_u = 0;
    this->plan_spin_mask_stable_p = 0;
    for (int stable_pos = 0; stable_pos < order; ++stable_pos) {
      int sp_u = this->plan_inv_argsort_u[stable_pos];
      int sp_p = this->plan_inv_argsort_p[stable_pos];
      if (this->master_unprimed.ops[sp_u].get_orbital_index() >= 1) this->plan_spin_mask_stable_u |= (1ULL << stable_pos);
      if (this->master_primed.ops[sp_p].get_orbital_index() >= 1) this->plan_spin_mask_stable_p |= (1ULL << stable_pos);
    }

    uint64_t full_mask = (order == 64) ? ~0ULL : (1ULL << order) - 1;
    plan.root_id       = this->solve_record(full_mask, full_mask, plan);
  }

  // ============================================================================
  //  Plan evaluation
  // ============================================================================

  namespace {
    std::vector<int> invert_argsort(std::vector<int> const &argsort) {
      std::vector<int> inv(argsort.size());
      for (size_t i = 0; i < argsort.size(); ++i) inv[argsort[i]] = (int)i;
      return inv;
    }

    template <typename T>
    Args<T> build_leaf_args_stable(CumulantPlan::LeafOps const &leaf, Args<T> const &master_unprimed, Args<T> const &master_primed,
                                   std::vector<int> const &inv_argsort_u, std::vector<int> const &inv_argsort_p) {
      int n = (int)leaf.u_global_idx.size();
      std::vector<double> taus;
      std::vector<FermionOperator<T>> ops;
      taus.reserve(2 * n);
      ops.reserve(2 * n);
      for (int i = 0; i < n; ++i) {
        int sp = inv_argsort_p[leaf.p_global_idx[i]];
        int su = inv_argsort_u[leaf.u_global_idx[i]];
        taus.push_back(master_primed.taus[sp]);
        ops.push_back(master_primed.ops[sp]);
        taus.push_back(master_unprimed.taus[su]);
        ops.push_back(master_unprimed.ops[su]);
      }
      return Args<T>(std::move(taus), std::move(ops));
    }
  } // namespace

  template <typename T>
  T evaluate_plan(CumulantPlan const &plan, Args<T> const &master_unprimed, Args<T> const &master_primed, HubbardSolver<T> const &solver,
                  bool infinite_U) {

    if (plan.root_id < 0) return T(0.0);

    std::vector<int> inv_argsort_u = invert_argsort(master_unprimed.argsort);
    std::vector<int> inv_argsort_p = invert_argsort(master_primed.argsort);

    std::vector<T> value;
    value.reserve(plan.nodes.size());

    for (size_t i = 0; i < plan.nodes.size(); ++i) {
      auto const &node = plan.nodes[i];

      Args<T> args = build_leaf_args_stable<T>(node.leaf, master_unprimed, master_primed, inv_argsort_u, inv_argsort_p);
      T v          = infinite_U ? solver.G0n_infinite_U(args) : solver.G0n(args);

      for (auto const &term : node.subtraction_terms) {
        T prod = T((double)term.sign);
        for (int fid : term.factor_node_ids) prod = prod * value[fid];
        v = v + prod;
      }
      value.push_back(v);
    }

    T C_stable = value[plan.root_id];
    return C_stable * T(master_unprimed.permutation_sign) * T(master_primed.permutation_sign);
  }

  // Explicit instantiations
  template class CumulantSolver<double>;
  template class CumulantSolver<Dual>;

  template double evaluate_plan<double>(CumulantPlan const &, Args<double> const &, Args<double> const &, HubbardSolver<double> const &, bool);
  template Dual evaluate_plan<Dual>(CumulantPlan const &, Args<Dual> const &, Args<Dual> const &, HubbardSolver<Dual> const &, bool);

} // namespace sc_expansion
