#include "vertex.hpp"
#include <unordered_map>

namespace sc_expansion {

  // Identify which (u_stable, p_stable) index pairs come from the same hopping
  // line (i.e. a self-loop: source and destination both at this vertex). The
  // cumulant plan uses this to keep the (c, c†) of a density insertion inside
  // the same partition block, so the vertex evaluates the density-density
  // cumulant instead of the 4-operator cumulant that would otherwise over-
  // subtract an "atomic propagator" off-diagonal term. Pre-sort positions in
  // the split unprimed/primed arrays coincide with stable indices, so walking
  // op_ids in input order is enough.
  template <typename T>
  static std::vector<std::pair<int, int>>
  compute_self_loop_pairs(std::vector<int> const &tau_indices, std::vector<uint8_t> const &op_ids) {
    std::unordered_map<int, int> line_count;
    for (int lidx : tau_indices) line_count[lidx]++;

    std::unordered_map<int, std::pair<int, int>> line_to_positions;
    int u_pos = 0;
    int p_pos = 0;
    for (size_t i = 0; i < op_ids.size(); ++i) {
      FermionOperator<T> f(op_ids[i]);
      auto &entry = line_to_positions[tau_indices[i]];
      if (f.get_action() == 0) {
        entry.first = u_pos++;
      } else {
        entry.second = p_pos++;
      }
    }

    std::vector<std::pair<int, int>> pairs;
    for (auto const &[lidx, cnt] : line_count) {
      if (cnt == 2) pairs.push_back(line_to_positions[lidx]);
    }
    return pairs;
  }

  template <typename T> VertexType<T>::VertexType(int n_legs_) : n_legs(n_legs_) {}

  template <typename T>
  VertexInstance<T>::VertexInstance(VertexType<T> *type_, std::vector<int> tau_indices_, std::vector<uint8_t> op_ids_)
     : type(type_), tau_indices(std::move(tau_indices_)), op_ids(std::move(op_ids_)) {}

  template <typename T>
  T VertexInstance<T>::get_value(const std::vector<double> &global_taus, const HubbardSolver<T> &solver, bool infinite_U) const {
    bool &dirty = infinite_U ? this->is_dirty_infinite : this->is_dirty_finite;
    T &cache    = infinite_U ? this->local_cache_infinite : this->local_cache_finite;

    if (!dirty) {
      this->type->record_local_hit();
      return cache;
    }
    this->type->record_local_miss();

    std::vector<double> local_taus;
    for (int idx : this->tau_indices) local_taus.push_back(global_taus[idx]);

    auto [unprimed, primed] = Args<T>::split_from_raw(local_taus, this->op_ids);

    if (!this->plan_built) {
      std::vector<double> dummy_taus(this->op_ids.size(), 0.5);
      auto [u0, p0] = Args<T>::split_from_raw(dummy_taus, this->op_ids);
      auto self_loop_pairs = compute_self_loop_pairs<T>(this->tau_indices, this->op_ids);
      CumulantSolver<T> builder(u0, p0, std::move(self_loop_pairs));
      builder.record_plan(this->plan);
      this->plan_built = true;
    }

    T canonical_val = evaluate_plan(this->plan, unprimed, primed, solver, infinite_U);

    cache = canonical_val * T(unprimed.permutation_sign) * T(primed.permutation_sign);
    dirty = false;

    return cache;
  }

  // Explicit instantiations
  template class VertexType<double>;
  template class VertexType<Dual>;

  template class VertexInstance<double>;
  template class VertexInstance<Dual>;

} // namespace sc_expansion
