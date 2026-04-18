#include "vertex.hpp"

namespace sc_expansion {

  template <int N_sites, typename T> VertexType<N_sites, T>::VertexType(int n_legs_) : n_legs(n_legs_) {}

  template <int N_sites, typename T>
  T VertexType<N_sites, T>::evaluate_canonical(Args<N_sites, T> const &unprimed, Args<N_sites, T> const &primed,
                                               HubbardSolver<N_sites, T> const &solver, bool infinite_U) const {

    // 1. Create the Canonical Key (The Args objects are ALREADY sorted by their constructor)
    std::vector<uint8_t> ops_u, ops_p;
    for (auto const &op : unprimed.ops) ops_u.push_back(op.op);
    for (auto const &op : primed.ops) ops_p.push_back(op.op);

    CanonicalVertexKey<N_sites, T> key{unprimed.taus, primed.taus, ops_u, ops_p, infinite_U};

    // 2. Check the Shared Spreadsheet
    auto it = this->global_cache.find(key);
    if (it != this->global_cache.end()) {
      // Hit! return value. (The Instance handles the Diagram-specific sign)
      this->cache_hits++;
      return it->second;
    }

    // 3. Cache Miss: Perform the expensive recursion
    this->cache_misses++;
    T result = compute_cumulant_decomposition(unprimed, primed, solver, infinite_U);

    // 4. Store for everyone else to use
    this->global_cache[key] = result;
    return result;
  }

  template <int N_sites, typename T>
  VertexInstance<N_sites, T>::VertexInstance(VertexType<N_sites, T> *type_, std::vector<int> tau_indices_, std::vector<uint8_t> op_ids_)
     : type(type_), tau_indices(std::move(tau_indices_)), op_ids(std::move(op_ids_)) {}

  template <int N_sites, typename T>
  T VertexInstance<N_sites, T>::get_value(const std::vector<double> &global_taus, const HubbardSolver<N_sites, T> &solver, bool infinite_U) const {
    // 1. Check the "Sticky Note" (Local Cache) — separate caches for finite/infinite-U
    bool &dirty = infinite_U ? this->is_dirty_infinite : this->is_dirty_finite;
    T &cache    = infinite_U ? this->local_cache_infinite : this->local_cache_finite;

    if (!dirty) return cache;

    // 2. Build local view from global taus
    std::vector<double> local_taus;
    for (int idx : this->tau_indices) local_taus.push_back(global_taus[idx]);

    // 3. Split and Sort (via Args)
    auto [unprimed, primed] = Args<N_sites, T>::split_from_raw(local_taus, this->op_ids);

    // 4. Build plans lazily on first evaluation (op pattern is fixed for the instance).
    if (!this->plans_built) {
      std::vector<double> dummy_taus(this->op_ids.size(), 0.5);
      auto [u0, p0]     = Args<N_sites, T>::split_from_raw(dummy_taus, this->op_ids);
      CumulantSolver<N_sites, T> builder_f(u0, p0, solver, /*infinite_U=*/false);
      builder_f.record_plan(this->plan_finite);
      CumulantSolver<N_sites, T> builder_i(u0, p0, solver, /*infinite_U=*/true);
      builder_i.record_plan(this->plan_infinite);
      this->plans_built = true;
    }

    // 5. Evaluate via the plan. Result is τ-sorted basis (same convention as the old
    // evaluate_canonical path), so the permutation_sign conversion below is unchanged.
    CumulantPlan const &plan = infinite_U ? this->plan_infinite : this->plan_finite;
    T canonical_val          = evaluate_plan(plan, unprimed, primed, solver, infinite_U);

    // 6. Restore the Diagram's leg convention using the saved signs
    // C(raw) = C(canonical) * sign(unprimed_sort) * sign(primed_sort)
    cache = canonical_val * T(unprimed.permutation_sign) * T(primed.permutation_sign);
    dirty = false;

    return cache;
  }

  // Explicit instantiations
  template class VertexType<1, double>;
  template class VertexType<1, Dual>;
  template class VertexType<2, double>;
  template class VertexType<2, Dual>;

  template class VertexInstance<1, double>;
  template class VertexInstance<1, Dual>;
  template class VertexInstance<2, double>;
  template class VertexInstance<2, Dual>;

} // namespace sc_expansion
