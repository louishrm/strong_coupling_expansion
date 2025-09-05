#include "./atomic_cumulant.hpp"

namespace atomic_cumulant {

  // A helper function to format a vector of any type into a comma-separated string.
  template <typename T> std::string format_vector(std::vector<T> const &vec) {
    std::stringstream ss;
    for (size_t i = 0; i < vec.size(); ++i) { ss << vec[i] << (i == vec.size() - 1 ? "" : ","); }
    return ss.str();
  }

  // Creates the string representation for a Green's function, e.g., "G(1,2 | 1',2')".
  std::string green_function(std::vector<std::string> const &unprimed, std::vector<std::string> const &primed) {
    std::stringstream ss;
    ss << "G(" << format_vector<std::string>(unprimed) << " | " << format_vector<std::string>(primed) << ")";
    return ss.str();
  }

  // Type definitions for clarity
  using Subset        = std::vector<std::string>;
  using Partition     = std::vector<Subset>;
  using AllPartitions = std::vector<Partition>;

  // Recursive helper function to generate all partitions of a set.
  void generate_partitions_recursive(const std::vector<std::string> &set, int index, Partition &current_partition, AllPartitions &result) {

    // Base case: If all elements have been placed into subsets, we have a complete partition.
    if (index == set.size()) {
      result.push_back(current_partition);
      return;
    }

    const std::string &element = set[index];

    // Choice 1: Add the current element to an existing subset in the current partition.
    for (int i = 0; i < current_partition.size(); i++) {
      current_partition[i].push_back(element);
      generate_partitions_recursive(set, index + 1, current_partition, result);
      current_partition[i].pop_back(); // Backtrack to explore other possibilities
    }

    // Choice 2: Create a new subset containing only the current element.
    current_partition.push_back({element});
    generate_partitions_recursive(set, index + 1, current_partition, result);
    current_partition.pop_back(); // Backtrack
  }

  // Main function to generate all partitions for a given set of strings.
  AllPartitions allPartitions(const std::vector<std::string> &set) {
    AllPartitions result;
    Partition current_partition; // This holds the partition being built during recursion
    generate_partitions_recursive(set, 0, current_partition, result);
    return result;
  }

  // Main function to perform and display the cumulant decomposition.
  void showCumulantDecomposition(std::vector<std::string> const &unprimed, std::vector<std::string> const &primed) {

    // Ensure the number of primed and unprimed indices are the same.
    if (unprimed.size() != primed.size() || unprimed.empty()) {
      std::cerr << "Error: Unprimed and primed sets must be non-empty and have the same size." << std::endl;
      return;
    }

    // The decomposition of a first-order cumulant is just the Green's function itself.
    if (unprimed.size() == 1) {
      std::cout << "C(" << unprimed[0] << "|" << primed[0] << ") = " << green_function(unprimed, primed) << std::endl;
      return;
    }

    // Generate all partitions of the unprimed set
    AllPartitions partitions = allPartitions(unprimed);

    // Create a mutable copy of the primed indices to generate permutations
    std::vector<std::string> current_primed_perm = primed;
    // Sort it first, as std::next_permutation generates permutations in lexicographical order
    std::sort(current_primed_perm.begin(), current_primed_perm.end());

    std::vector<std::string> final_terms;

    // --- Core Logic ---
    // Loop over all unique permutations of the primed indices.
    do {
      // For each permutation, loop over all partitions of the unprimed indices.
      for (auto const &partition : partitions) {

        std::stringstream term_ss;
        int primed_cursor = 0; // This acts as a cursor for the current_primed_perm vector.

        // For each subset in the current partition...
        for (size_t i = 0; i < partition.size(); ++i) {
          auto const &subset = partition[i];

          // Create the matching subset of primed indices by taking a "slice"
          // of the current permutation. The size of the slice is determined by the unprimed subset's size.
          std::vector<std::string> primed_match(current_primed_perm.begin() + primed_cursor,
                                                current_primed_perm.begin() + primed_cursor + subset.size());

          // As per the convention, the primed indices within a Green's function should be ordered.
          std::sort(primed_match.begin(), primed_match.end());

          // Advance the cursor for the next subset.
          primed_cursor += subset.size();

          // Append the corresponding Green's function to the term.
          term_ss << green_function(subset, primed_match);
          if (i < partition.size() - 1) {
            term_ss << " * "; // Products of GFs form a single term in the sum
          }
        }
        final_terms.push_back(term_ss.str());
      }
    } while (std::next_permutation(current_primed_perm.begin(), current_primed_perm.end()));

    // Enforcing a canonical order on primed indices can lead to duplicate terms.
    // For example, G(1,2|2',1') becomes G(1,2|1',2'), which might be generated by another permutation.
    // We sort all generated term strings and remove duplicates to get the final unique set.
    std::sort(final_terms.begin(), final_terms.end());
    final_terms.erase(std::unique(final_terms.begin(), final_terms.end()), final_terms.end());

    // --- Final Output ---
    // Print the final, formatted decomposition.
    std::cout << "C(" << format_vector(unprimed) << " | " << format_vector(primed) << ") = \n";
    for (size_t i = 0; i < final_terms.size(); ++i) { std::cout << "  " << final_terms[i] << (i == final_terms.size() - 1 ? "" : " + \n"); }
    std::cout << std::endl;
  }

} // namespace atomic_cumulant

// Template instantiations (if you compile this as a separate .cpp file)
template std::string atomic_cumulant::format_vector<std::string>(std::vector<std::string> const &);
template std::string atomic_cumulant::format_vector<int>(std::vector<int> const &);