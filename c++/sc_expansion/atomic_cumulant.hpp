#ifndef ATOMIC_CUMULANT_HPP
#define ATOMIC_CUMULANT_HPP

#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

namespace atomic_cumulant {

  // --- Type Aliases ---
  using Subset        = std::vector<std::string>;
  using Partition     = std::vector<Subset>;
  using AllPartitions = std::vector<Partition>;

  // --- Template Function Definition ---
  // The full definition of a template function must be in the header file.
  // Formats a vector of any type into a comma-separated string.
  template <typename T> std::string format_vector(std::vector<T> const &vec);

  // Formats a Green's function string representation.
  std::string green_function(std::vector<std::string> const &unprimed, std::vector<std::string> const &primed);

  // Generates all possible partitions of a given set of strings.
  AllPartitions allPartitions(const std::vector<std::string> &set);

  void showCumulantDecomposition(std::vector<std::string> const &unprimed, std::vector<std::string> const &primed);

} // namespace atomic_cumulant

#endif // ATOMIC_CUMULANT_HPP