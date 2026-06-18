#include "cumulant.hpp"
#include "dual.hpp"
#include <numeric>
#include <algorithm>
#include <functional>
#include <iostream>

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

  template <int N_sites, typename T>
  CumulantSolver<N_sites, T>::CumulantSolver(Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed)
     : master_unprimed(unprimed), master_primed(primed) {}

  template <int N_sites, typename T>
  CumulantSolver<N_sites, T>::CumulantSolver(Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed,
                                              HubbardSolver<N_sites, T> const &solver, bool infinite_U_)
     : master_unprimed(unprimed), master_primed(primed), solver_ptr(&solver), infinite_U(infinite_U_) {}

  template <int N_sites, typename T> T CumulantSolver<N_sites, T>::compute_cumulant_decomposition() {
    CumulantPlan plan;
    this->record_plan(plan);
    return evaluate_plan(plan, this->master_unprimed, this->master_primed, *this->solver_ptr, this->infinite_U);
  }

  template <int N_sites, typename T>
  void CumulantSolver<N_sites, T>::record_distribute_primed_collect(
     std::vector<uint64_t> const &u_partition_masks, int u_idx, uint64_t current_p_pool, int overall_sign,
     CumulantSolver::OpBlockList &block_list_so_far, std::vector<int> const &stable_map_u, std::vector<int> const &stable_map_p,
     std::vector<std::pair<int, CumulantSolver::OpBlockList>> &out,
     std::vector<int> const &special_u_idxs, std::vector<int> const &local_p_block_idxs) {

    if (u_idx == (int)u_partition_masks.size()) {
      out.emplace_back(overall_sign, block_list_so_far);
      return;
    }

    uint64_t u_mask = u_partition_masks[u_idx];
    int needed_k    = popcount(u_mask);

    for_each_subset(current_p_pool, needed_k, [&](uint64_t p_submask) {
      // Block-constraint filter: for each active block k, the U-subset that
      // contains the mark unprimed leg (special_u_idxs[k]) MUST be paired
      // with a P-subset containing the mark primed leg (local_p_block_idxs[k]);
      // every other U-subset MUST exclude that primed leg. Different blocks
      // act independently — they may share or split subsets freely.
      for (size_t k = 0; k < special_u_idxs.size(); ++k) {
        bool p_submask_has_block = (p_submask >> local_p_block_idxs[k]) & 1ULL;
        bool must_have           = (u_idx == special_u_idxs[k]);
        if (must_have != p_submask_has_block) return;
      }
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

      block_list_so_far.emplace_back(global_mask_u_stable, global_mask_p_stable);
      this->record_distribute_primed_collect(u_partition_masks, u_idx + 1, current_p_pool ^ p_submask, overall_sign * step_sign_p, block_list_so_far,
                                             stable_map_u, stable_map_p, out, special_u_idxs, local_p_block_idxs);
      block_list_so_far.pop_back();
    });
  }

  // Compact 4-bits-per-orbital count encoding for the decoration multiset.
  // 0 (empty) preserves the original cache-key shape for un-decorated nodes.
  static uint64_t encode_decoration(std::vector<int> const &decoration) {
    uint64_t k = 0;
    for (int sigma : decoration) { k += (1ULL << (4 * sigma)); }
    return k;
  }

  template <int N_sites, typename T>
  void CumulantSolver<N_sites, T>::enumerate_atom_routings_and_emit(
     std::vector<int> const &decoration, int atom_idx, CumulantSolver::OpBlockList const &op_blocks, int base_sign,
     std::vector<std::vector<int>> &per_op_block_decoration, std::vector<std::vector<int>> &free_groups,
     std::vector<CumulantPlan::ProductTerm> &out_terms, CumulantPlan &plan) {

    if (atom_idx == (int)decoration.size()) {
      // Skip the trivial joint partition (single block = the leaf itself).
      int total_blocks = (int)op_blocks.size() + (int)free_groups.size();
      if (total_blocks <= 1) return;

      CumulantPlan::ProductTerm term;
      term.sign = base_sign;
      for (size_t k = 0; k < op_blocks.size(); ++k) {
        int child = this->solve_record(op_blocks[k].first, op_blocks[k].second, per_op_block_decoration[k], plan);
        if (child < 0) return;
        term.factor_node_ids.push_back(child);
      }
      for (auto const &fg : free_groups) {
        int child = this->solve_record(/*mask_u=*/0, /*mask_p=*/0, fg, plan);
        if (child < 0) return;
        term.factor_node_ids.push_back(child);
      }
      out_terms.push_back(std::move(term));
      return;
    }

    int sigma = decoration[atom_idx];
    // Option 1: attach this atom to op-block k.
    for (size_t k = 0; k < op_blocks.size(); ++k) {
      per_op_block_decoration[k].push_back(sigma);
      this->enumerate_atom_routings_and_emit(decoration, atom_idx + 1, op_blocks, base_sign, per_op_block_decoration, free_groups, out_terms, plan);
      per_op_block_decoration[k].pop_back();
    }
    // Option 2: attach to an existing free group g.
    for (size_t g = 0; g < free_groups.size(); ++g) {
      free_groups[g].push_back(sigma);
      this->enumerate_atom_routings_and_emit(decoration, atom_idx + 1, op_blocks, base_sign, per_op_block_decoration, free_groups, out_terms, plan);
      free_groups[g].pop_back();
    }
    // Option 3: open a new free group containing this atom alone.
    free_groups.push_back({sigma});
    this->enumerate_atom_routings_and_emit(decoration, atom_idx + 1, op_blocks, base_sign, per_op_block_decoration, free_groups, out_terms, plan);
    free_groups.pop_back();
  }

  template <int N_sites, typename T>
  int CumulantSolver<N_sites, T>::solve_record(uint64_t mask_u_stable, uint64_t mask_p_stable, std::vector<int> const &decoration,
                                                CumulantPlan &plan) {

    if (popcount(mask_u_stable & this->plan_spin_mask_stable_u) != popcount(mask_p_stable & this->plan_spin_mask_stable_p)) return -1;

    CacheKey key{mask_u_stable, mask_p_stable, encode_decoration(decoration)};
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
    node.leaf.u_global_idx     = stable_map_u;
    node.leaf.p_global_idx     = stable_map_p;
    node.leaf_density_orbitals = decoration;

    // Base cases (no further decomposition possible):
    //   order == 1 + no decoration: existing rooted-mark size-1 sub-cumulant.
    //   order == 0 + |decoration| <= 1: leaf-only (no atoms = empty trace, |1 atom| = ⟨n_σ⟩).
    // Same-spin idempotency (n_σ^k = n_σ) is handled by G0n_with_densities at eval time.
    if ((order == 1 && decoration.empty()) || (order == 0 && decoration.size() <= 1)) {
      int new_id = (int)plan.nodes.size();
      plan.nodes.push_back(std::move(node));
      this->plan_node_ids.emplace(key, new_id);
      return new_id;
    }

    // Per-block constraint state (existing add_block_constraint mechanism).
    // Independent of, and orthogonal to, the new density-decoration mechanism.
    std::vector<int> active_local_u_block_idxs;
    std::vector<int> active_local_p_block_idxs;
    active_local_u_block_idxs.reserve(this->block_u_stable_list.size());
    active_local_p_block_idxs.reserve(this->block_p_stable_list.size());
    for (size_t k = 0; k < this->block_u_stable_list.size(); ++k) {
      int bu = this->block_u_stable_list[k];
      int bp = this->block_p_stable_list[k];
      bool present = (mask_u_stable & (1ULL << bu)) && (mask_p_stable & (1ULL << bp));
      if (!present) continue;
      active_local_u_block_idxs.push_back(popcount(mask_u_stable & ((1ULL << bu) - 1)));
      active_local_p_block_idxs.push_back(popcount(mask_p_stable & ((1ULL << bp) - 1)));
    }

    auto const &unprimed_partitions = get_partitions(order);
    for (auto const &partition : unprimed_partitions) {
      // Skip the trivial u-partition only when no atoms can split off; with a
      // non-empty decoration we still need to enumerate atom routings even on
      // the trivial op-partition (atoms can break off into free groups).
      if (partition.size() == 1 && decoration.empty()) continue;

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

      std::vector<int> special_u_idxs(active_local_u_block_idxs.size(), -1);
      for (size_t k = 0; k < active_local_u_block_idxs.size(); ++k) {
        int lub = active_local_u_block_idxs[k];
        for (size_t s = 0; s < u_partition_local_masks.size(); ++s) {
          if (u_partition_local_masks[s] & (1ULL << lub)) { special_u_idxs[k] = (int)s; break; }
        }
      }

      int base_sign              = -sign_u;
      uint64_t full_local_p_mask = (order == 64) ? ~0ULL : (1ULL << order) - 1;

      // Collect all (overall_sign, op-block-list) tuples for this u-partition.
      std::vector<std::pair<int, OpBlockList>> collected;
      OpBlockList block_list_scratch;
      block_list_scratch.reserve(partition.size());
      this->record_distribute_primed_collect(u_partition_local_masks, 0, full_local_p_mask, base_sign, block_list_scratch, stable_map_u, stable_map_p,
                                             collected, special_u_idxs, active_local_p_block_idxs);

      // For each op-partition, enumerate per-atom routings on top.
      for (auto const &entry : collected) {
        int op_sign                       = entry.first;
        OpBlockList const &op_blocks      = entry.second;
        std::vector<std::vector<int>> per_op_block_decoration(op_blocks.size());
        std::vector<std::vector<int>> free_groups;
        this->enumerate_atom_routings_and_emit(decoration, /*atom_idx=*/0, op_blocks, op_sign, per_op_block_decoration, free_groups,
                                               node.subtraction_terms, plan);
      }
    }

    int new_id = (int)plan.nodes.size();
    plan.nodes.push_back(std::move(node));
    this->plan_node_ids.emplace(key, new_id);
    return new_id;
  }

  template <int N_sites, typename T> void CumulantSolver<N_sites, T>::record_plan(CumulantPlan &plan) {
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
      if (this->master_unprimed.ops[sp_u].get_orbital_index() >= N_sites) this->plan_spin_mask_stable_u |= (1ULL << stable_pos);
      if (this->master_primed.ops[sp_p].get_orbital_index() >= N_sites) this->plan_spin_mask_stable_p |= (1ULL << stable_pos);
    }

    uint64_t full_mask = (order == 64) ? ~0ULL : (1ULL << order) - 1;
    plan.root_id       = this->solve_record(full_mask, full_mask, this->pending_static_densities, plan);

    // For each block constraint, locate the size-1 sub-cumulant node whose
    // only operators are that block's mark legs. Its leaf G0n would be
    // ⟨T_τ c_σ(0) c†_σ(0)⟩ — convention-dependent at equal times — so we
    // override it with the unambiguous ⟨n_σ⟩ from the solver at eval time.
    plan.mark_node_ids.clear();
    plan.mark_orbitals.clear();
    for (size_t k = 0; k < this->block_u_stable_list.size(); ++k) {
      int bu = this->block_u_stable_list[k];
      int bp = this->block_p_stable_list[k];
      CacheKey mark_key{1ULL << bu, 1ULL << bp, /*decoration_key=*/0};
      if (auto it = this->plan_node_ids.find(mark_key); it != this->plan_node_ids.end()) {
        plan.mark_node_ids.push_back(it->second);
        int sp_u = this->plan_inv_argsort_u[bu];
        plan.mark_orbitals.push_back((int)this->master_unprimed.ops[sp_u].get_orbital_index());
      }
    }

    // Record coincidence groups, expanding block indices to per-leg stable
    // indices. Then scan all non-trivial subsets (size >= 2) of each group
    // and register any plan node whose mask is exactly the subset's union
    // as a mark override returning ⟨n_σ⟩ — valid because n_σ(0)^k = n_σ(0).
    plan.coincidence_groups.clear();
    for (auto const &pg : this->pending_coincidence_groups) {
      CumulantPlan::CoincidenceGroup g;
      g.orbital = pg.orbital;
      for (int blk : pg.block_indices) {
        g.u_stable.push_back(this->block_u_stable_list[blk]);
        g.p_stable.push_back(this->block_p_stable_list[blk]);
      }
      plan.coincidence_groups.push_back(std::move(g));
    }

    for (auto const &g : plan.coincidence_groups) {
      int n_blk = (int)g.u_stable.size();
      if (n_blk < 2) continue;
      for (uint32_t sub = 1; sub < (1U << n_blk); ++sub) {
        if (__builtin_popcount(sub) < 2) continue;
        uint64_t u_mask = 0, p_mask = 0;
        for (int j = 0; j < n_blk; ++j) {
          if ((sub >> j) & 1U) {
            u_mask |= 1ULL << g.u_stable[j];
            p_mask |= 1ULL << g.p_stable[j];
          }
        }
        CacheKey key{u_mask, p_mask, /*decoration_key=*/0};
        if (auto it = this->plan_node_ids.find(key); it != this->plan_node_ids.end()) {
          plan.mark_node_ids.push_back(it->second);
          plan.mark_orbitals.push_back(g.orbital);
        }
      }
    }
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

    template <int N_sites, typename T>
    Args<N_sites, T> build_leaf_args_stable(CumulantPlan::LeafOps const &leaf, Args<N_sites, T> const &master_unprimed,
                                            Args<N_sites, T> const &master_primed, std::vector<int> const &inv_argsort_u,
                                            std::vector<int> const &inv_argsort_p,
                                            std::vector<CumulantPlan::CoincidenceGroup> const &groups) {
      // For each coincidence group, identify which of its u/p-stable
      // indices are present in this leaf. If >= 2 marks of one group are
      // present, the product c†_σ(0)·c_σ(0)·… reduces to a single n_σ(0).
      // Drop all but the first occurrence in each list independently.
      // (Lists drop the same *count* since the block constraint pairs each
      // mark's u-leg with its p-leg in the same partition cell; positions
      // need not align — the leaf is a joint expectation of an unordered
      // operator set, not a position-paired list.)
      std::vector<bool> drop_u(leaf.u_global_idx.size(), false);
      std::vector<bool> drop_p(leaf.p_global_idx.size(), false);
      for (auto const &g : groups) {
        std::vector<int> present_u_pos, present_p_pos;
        for (size_t i = 0; i < leaf.u_global_idx.size(); ++i) {
          if (std::find(g.u_stable.begin(), g.u_stable.end(), leaf.u_global_idx[i]) != g.u_stable.end())
            present_u_pos.push_back((int)i);
        }
        for (size_t i = 0; i < leaf.p_global_idx.size(); ++i) {
          if (std::find(g.p_stable.begin(), g.p_stable.end(), leaf.p_global_idx[i]) != g.p_stable.end())
            present_p_pos.push_back((int)i);
        }
        if (present_u_pos.size() < 2) continue;
        for (size_t i = 1; i < present_u_pos.size(); ++i) drop_u[present_u_pos[i]] = true;
        for (size_t i = 1; i < present_p_pos.size(); ++i) drop_p[present_p_pos[i]] = true;
      }

      // Collect kept u and p separately, then re-pair (Args sorts by tau).
      std::vector<int> kept_u_global, kept_p_global;
      for (size_t i = 0; i < leaf.u_global_idx.size(); ++i)
        if (!drop_u[i]) kept_u_global.push_back(leaf.u_global_idx[i]);
      for (size_t i = 0; i < leaf.p_global_idx.size(); ++i)
        if (!drop_p[i]) kept_p_global.push_back(leaf.p_global_idx[i]);

      int n_eff = (int)kept_u_global.size();
      std::vector<double> taus;
      std::vector<FermionOperator<N_sites, T>> ops;
      taus.reserve(2 * n_eff);
      ops.reserve(2 * n_eff);
      for (int i = 0; i < n_eff; ++i) {
        int sp = inv_argsort_p[kept_p_global[i]];
        int su = inv_argsort_u[kept_u_global[i]];
        taus.push_back(master_primed.taus[sp]);
        ops.push_back(master_primed.ops[sp]);
        taus.push_back(master_unprimed.taus[su]);
        ops.push_back(master_unprimed.ops[su]);
      }
      return Args<N_sites, T>(std::move(taus), std::move(ops));
    }
  } // namespace

  template <int N_sites, typename T>
  T evaluate_plan(CumulantPlan const &plan, Args<N_sites, T> const &master_unprimed, Args<N_sites, T> const &master_primed,
                  HubbardSolver<N_sites, T> const &solver, bool infinite_U) {

    if (plan.root_id < 0) return T(0.0);

    std::vector<int> inv_argsort_u = invert_argsort(master_unprimed.argsort);
    std::vector<int> inv_argsort_p = invert_argsort(master_primed.argsort);

    std::vector<T> value;
    value.reserve(plan.nodes.size());

    for (size_t i = 0; i < plan.nodes.size(); ++i) {
      auto const &node = plan.nodes[i];

      T v;
      int mark_slot = -1;
      for (size_t k = 0; k < plan.mark_node_ids.size(); ++k) {
        if ((int)i == plan.mark_node_ids[k]) { mark_slot = (int)k; break; }
      }
      if (mark_slot >= 0) {
        // Mark sub-block at finite U: use ⟨n_σ⟩ directly rather than the
        // T_τ-ordered equal-time evaluation of (c_σ(0), c†_σ(0)). At
        // infinite U the closed-form `G0n_infinite_U` already returns
        // ⟨n_σ⟩_∞ for this equal-time pair, so fall through there.
        v = solver.compute_n_sigma(plan.mark_orbitals[mark_slot]);
      } else {
        Args<N_sites, T> args =
           build_leaf_args_stable<N_sites, T>(node.leaf, master_unprimed, master_primed, inv_argsort_u, inv_argsort_p, plan.coincidence_groups);
        if (node.leaf_density_orbitals.empty()) {
          v = infinite_U ? solver.G0n_infinite_U(args) : solver.G0n(args);
        } else {
          v = infinite_U ? solver.G0n_with_densities_infinite_U(args, node.leaf_density_orbitals)
                         : solver.G0n_with_densities(args, node.leaf_density_orbitals);
        }
      }

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

  template <int N_sites, typename T>
  T evaluate_plan_incremental(CumulantPlan const &plan, Args<N_sites, T> const &master_unprimed, Args<N_sites, T> const &master_primed,
                              HubbardSolver<N_sites, T> const &solver, bool infinite_U, std::vector<T> &value,
                              std::vector<uint64_t> const &node_line_mask, bool recompute_all, uint64_t changed_lines) {

    if (plan.root_id < 0) return T(0.0);

    value.resize(plan.nodes.size()); // no-op after the first call; guarantees sizing

    std::vector<int> inv_argsort_u = invert_argsort(master_unprimed.argsort);
    std::vector<int> inv_argsort_p = invert_argsort(master_primed.argsort);

    for (size_t i = 0; i < plan.nodes.size(); ++i) {
      // Reuse the cached value when neither this node's leaf nor any descendant
      // sub-cumulant depends on a changed line. value[i] holds the raw node
      // value (the master permutation sign is applied once at the return), so it
      // is reusable regardless of how the master sort moved for OTHER lines.
      if (!recompute_all && (node_line_mask[i] & changed_lines) == 0) continue;

      auto const &node = plan.nodes[i];

      T v;
      int mark_slot = -1;
      for (size_t k = 0; k < plan.mark_node_ids.size(); ++k) {
        if ((int)i == plan.mark_node_ids[k]) { mark_slot = (int)k; break; }
      }
      if (mark_slot >= 0) {
        v = solver.compute_n_sigma(plan.mark_orbitals[mark_slot]);
      } else {
        Args<N_sites, T> args =
           build_leaf_args_stable<N_sites, T>(node.leaf, master_unprimed, master_primed, inv_argsort_u, inv_argsort_p, plan.coincidence_groups);
        if (node.leaf_density_orbitals.empty()) {
          v = infinite_U ? solver.G0n_infinite_U(args) : solver.G0n(args);
        } else {
          v = infinite_U ? solver.G0n_with_densities_infinite_U(args, node.leaf_density_orbitals)
                         : solver.G0n_with_densities(args, node.leaf_density_orbitals);
        }
      }

      for (auto const &term : node.subtraction_terms) {
        T prod = T((double)term.sign);
        for (int fid : term.factor_node_ids) prod = prod * value[fid];
        v = v + prod;
      }
      value[i] = v;
    }

    T C_stable = value[plan.root_id];
    return C_stable * T(master_unprimed.permutation_sign) * T(master_primed.permutation_sign);
  }

  // Explicit instantiations
  template class CumulantSolver<1, double>;
  template class CumulantSolver<1, Dual>;
  template class CumulantSolver<2, double>;
  template class CumulantSolver<2, Dual>;

  template double evaluate_plan<1, double>(CumulantPlan const &, Args<1, double> const &, Args<1, double> const &,
                                           HubbardSolver<1, double> const &, bool);
  template Dual evaluate_plan<1, Dual>(CumulantPlan const &, Args<1, Dual> const &, Args<1, Dual> const &,
                                       HubbardSolver<1, Dual> const &, bool);
  template double evaluate_plan<2, double>(CumulantPlan const &, Args<2, double> const &, Args<2, double> const &,
                                           HubbardSolver<2, double> const &, bool);
  template Dual evaluate_plan<2, Dual>(CumulantPlan const &, Args<2, Dual> const &, Args<2, Dual> const &,
                                       HubbardSolver<2, Dual> const &, bool);

  template double evaluate_plan_incremental<1, double>(CumulantPlan const &, Args<1, double> const &, Args<1, double> const &,
                                                       HubbardSolver<1, double> const &, bool, std::vector<double> &, std::vector<uint64_t> const &,
                                                       bool, uint64_t);
  template Dual evaluate_plan_incremental<1, Dual>(CumulantPlan const &, Args<1, Dual> const &, Args<1, Dual> const &, HubbardSolver<1, Dual> const &,
                                                   bool, std::vector<Dual> &, std::vector<uint64_t> const &, bool, uint64_t);
  template double evaluate_plan_incremental<2, double>(CumulantPlan const &, Args<2, double> const &, Args<2, double> const &,
                                                       HubbardSolver<2, double> const &, bool, std::vector<double> &, std::vector<uint64_t> const &,
                                                       bool, uint64_t);
  template Dual evaluate_plan_incremental<2, Dual>(CumulantPlan const &, Args<2, Dual> const &, Args<2, Dual> const &, HubbardSolver<2, Dual> const &,
                                                   bool, std::vector<Dual> &, std::vector<uint64_t> const &, bool, uint64_t);

} // namespace sc_expansion
