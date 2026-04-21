#pragma once

#include <vector>
#include "args.hpp"
#include "cumulant.hpp"
#include "dual.hpp"

namespace sc_expansion {

  // VertexType acts as the "Instruction Manual" shared by all VertexInstances of a given leg count.
  template <typename T> class VertexType {
    public:
    explicit VertexType(int n_legs);

    std::pair<long, long> get_local_cache_stats() const { return {this->local_cache_hits, this->local_cache_misses}; }

    void record_local_hit() const { ++this->local_cache_hits; }
    void record_local_miss() const { ++this->local_cache_misses; }

    private:
    int n_legs;
    mutable long local_cache_hits   = 0;
    mutable long local_cache_misses = 0;
  };

  // VertexInstance is the specific "LEGO brick" in a Diagram
  template <typename T> class VertexInstance {
    public:
    VertexInstance(VertexType<T> *type_, std::vector<int> tau_indices_, std::vector<uint8_t> op_ids_);

    T get_value(const std::vector<double> &global_taus, const HubbardSolver<T> &solver, bool infinite_U) const;

    void mark_dirty() {
      this->is_dirty_finite   = true;
      this->is_dirty_infinite = true;
    }

    const std::vector<int> &get_tau_indices() const { return this->tau_indices; }

    private:
    VertexType<T> *type;
    std::vector<int> tau_indices;
    std::vector<uint8_t> op_ids;

    mutable T local_cache_finite;
    mutable T local_cache_infinite;
    mutable bool is_dirty_finite   = true;
    mutable bool is_dirty_infinite = true;

    mutable CumulantPlan plan;
    mutable bool plan_built = false;
  };

} // namespace sc_expansion
