#pragma once

#include <vector>
#include "args.hpp"
#include "cumulant.hpp"
#include "dual.hpp"

namespace sc_expansion {

  struct VertexCacheKey {
    std::vector<int> op_indices; // Specific actions (Site, Spin, Dag)
    std::vector<double> taus;    // Specific times from the MCMC
    bool infinite_U;             // Whether we are in the high-U limit

    bool operator==(const VertexCacheKey &other) const {
      return op_indices == other.op_indices && taus == other.taus && infinite_U == other.infinite_U; //times are always exactly preserved if unchanged
    }
  };

  struct VertexCacheHasher {
    std::size_t operator()(const VertexCacheKey &k) const {
      // Use a robust hash combine (e.g., boost::hash_combine style)
      std::size_t h = 0;
      for (int op : k.op_indices) h ^= std::hash<int>{}(op) + 0x9e3779b9 + (h << 6) + (h >> 2);
      for (double t : k.taus) h ^= std::hash<double>{}(t) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<bool>{}(k.infinite_U) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  template <int N_sites, typename T> class VertexType {

    public:
    VertexType(int n_legs);

    T evaluate(std::vector<double> const &taus, std::vector<int> const &op_indices, HubbardSolver<N_sites, T> const &solver, bool infinite_U) const;

    private:
    int n_legs;

    mutable std::unordered_map<VertexCacheKey, T, VertexCacheHasher> global_cache;

    std::pair<ArgList, ArgList> map_args_to_list(std::vector<double> const &taus, std::vector<int> const &op_indices) const;
  };

  class VertexInstance {};
} // namespace sc_expansion