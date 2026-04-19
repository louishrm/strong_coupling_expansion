#pragma once

#include <vector>
#include <unordered_map>
#include "args.hpp"
#include "cumulant.hpp"
#include "dual.hpp"

namespace sc_expansion {

  // The key for global memoization. Uses the CANONICAL (time-sorted) operator and time arrays.
  template <typename T> struct CanonicalVertexKey {
    std::vector<double> taus_u;
    std::vector<double> taus_p;
    std::vector<uint8_t> ops_u;
    std::vector<uint8_t> ops_p;
    bool infinite_U;

    bool operator==(const CanonicalVertexKey &other) const {
      return infinite_U == other.infinite_U && taus_u == other.taus_u && taus_p == other.taus_p && ops_u == other.ops_u && ops_p == other.ops_p;
    }
  };

  template <typename T> struct CanonicalVertexHasher {
    std::size_t operator()(const CanonicalVertexKey<T> &k) const {
      std::size_t h = 0;
      auto combine  = [&](auto const &v) {
        for (auto const &x : v) h ^= std::hash<std::decay_t<decltype(x)>>{}(x) + 0x9e3779b9 + (h << 6) + (h >> 2);
      };
      combine(k.taus_u);
      combine(k.taus_p);
      combine(k.ops_u);
      combine(k.ops_p);
      h ^= std::hash<bool>{}(k.infinite_U) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  // VertexType acts as the "Instruction Manual" and the "Shared Spreadsheet"
  template <typename T> class VertexType {
    public:
    explicit VertexType(int n_legs);

    std::pair<long, long> get_cache_stats() const { return {this->cache_hits, this->cache_misses}; }
    size_t get_cache_size() const { return this->global_cache.size(); }
    void clear_global_cache() const { this->global_cache.clear(); }

    private:
    int n_legs;
    mutable std::unordered_map<CanonicalVertexKey<T>, T, CanonicalVertexHasher<T>> global_cache;
    mutable long cache_hits   = 0;
    mutable long cache_misses = 0;
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
