#pragma once
#include <vector>
#include <algorithm>
#include <numeric>

namespace sc_expansion {

  inline uint64_t factorial(uint64_t k) {
    if (k <= 1) return 1;
    return k * factorial(k - 1);
  }

  inline std::vector<std::vector<int>> generate_permutations(int n) {

    std::vector<std::vector<int>> permutations;
    std::vector<int> vertices(n);
    for (int i = 0; i < n; ++i) vertices[i] = i;

    std::sort(vertices.begin(), vertices.end());
    do { permutations.push_back(vertices); } while (std::next_permutation(vertices.begin(), vertices.end()));
    return permutations;
  }

} // namespace sc_expansion