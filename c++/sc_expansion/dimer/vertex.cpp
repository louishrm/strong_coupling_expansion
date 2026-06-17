#include "vertex.hpp"
#include "../dual.hpp"

namespace sc_expansion::dimer {

  template <typename T> VertexType<T>::VertexType(int n_legs_) : n_legs(n_legs_) {}

  template <typename T>
  VertexInstance<T>::VertexInstance(VertexType<T> *type_, std::vector<int> tau_indices_, std::vector<uint8_t> op_ids_)
     : type(type_), tau_indices(std::move(tau_indices_)), op_ids(std::move(op_ids_)) {}

  template <typename T>
  T VertexInstance<T>::get_value(std::vector<double> const &global_taus, HubbardSolver<2, T> const &solver) const {
    if (!this->is_dirty_finite) {
      this->type->record_local_hit();
      return this->local_cache_finite;
    }
    this->type->record_local_miss();

    std::vector<double> local_taus;
    for (int idx : this->tau_indices) local_taus.push_back(global_taus[idx]);

    auto [unprimed, primed] = Args<2, T>::split_from_raw(local_taus, this->op_ids);

    if (!this->plan_built) {
      std::vector<double> dummy_taus(this->op_ids.size(), 0.5);
      auto [u0, p0] = Args<2, T>::split_from_raw(dummy_taus, this->op_ids);
      CumulantSolver<2, T> builder(u0, p0);
      builder.record_plan(this->plan);
      this->plan_built = true;
    }

    T canonical_val = evaluate_plan(this->plan, unprimed, primed, solver, /*infinite_U=*/false);

    this->local_cache_finite = canonical_val * T(unprimed.permutation_sign) * T(primed.permutation_sign);
    this->is_dirty_finite    = false;

    return this->local_cache_finite;
  }

  template class VertexType<double>;
  template class VertexType<Dual>;

  template class VertexInstance<double>;
  template class VertexInstance<Dual>;

} // namespace sc_expansion::dimer
