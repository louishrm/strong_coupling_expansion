#include "vertex.hpp"
#include "../dual.hpp"

namespace sc_expansion::atomic {

  template <typename T> VertexType<T>::VertexType(int n_legs_) : n_legs(n_legs_) {}

  template <typename T>
  VertexInstance<T>::VertexInstance(VertexType<T> *type_, std::vector<int> tau_indices_, std::vector<uint8_t> op_ids_)
     : type(type_), tau_indices(std::move(tau_indices_)), op_ids(std::move(op_ids_)) {}

  template <typename T>
  T VertexInstance<T>::get_value(std::vector<double> const &global_taus, HubbardSolver<1, T> const &solver, bool infinite_U) const {
    bool &dirty = infinite_U ? this->is_dirty_infinite : this->is_dirty_finite;
    T &cache    = infinite_U ? this->local_cache_infinite : this->local_cache_finite;

    if (!dirty) {
      this->type->record_local_hit();
      return cache;
    }
    this->type->record_local_miss();

    std::vector<double> local_taus;
    for (int idx : this->tau_indices) local_taus.push_back(global_taus[idx]);

    auto [unprimed, primed] = Args<1, T>::split_from_raw(local_taus, this->op_ids);

    if (!this->plan_built) {
      std::vector<double> dummy_taus(this->op_ids.size(), 0.5);
      auto [u0, p0] = Args<1, T>::split_from_raw(dummy_taus, this->op_ids);
      CumulantSolver<1, T> builder(u0, p0);
      for (size_t k = 0; k < this->block_u_inputs.size(); ++k) {
        builder.add_block_constraint(this->block_u_inputs[k], this->block_p_inputs[k]);
      }
      for (size_t k = 0; k < this->coincidence_orbitals.size(); ++k) {
        builder.add_coincidence_group(this->coincidence_orbitals[k], this->coincidence_block_indices[k]);
      }
      for (int orbital : this->static_density_orbitals) { builder.add_static_density(orbital); }
      builder.record_plan(this->plan);
      this->plan_built = true;
    }

    T canonical_val = evaluate_plan(this->plan, unprimed, primed, solver, infinite_U);

    cache = canonical_val * T(unprimed.permutation_sign) * T(primed.permutation_sign);
    dirty = false;

    return cache;
  }

  template class VertexType<double>;
  template class VertexType<Dual>;

  template class VertexInstance<double>;
  template class VertexInstance<Dual>;

} // namespace sc_expansion::atomic
