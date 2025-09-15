#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <sstream>

// Type definitions for clarity
using subset         = std::vector<int>;
using partition      = std::vector<subset>;
using all_partitions = std::vector<partition>;

void get_partitions(std::vector<int> &set, int index, all_partitions &ans, partition &current_partition) {
  //If all elements are processed, then finished, add the current partition.
  if (index == set.size()) {
    ans.push_back(current_partition);
    return;
  }

  // For each subset in the partition
  // add the current element to it
  // and recall
  for (int i = 0; i < current_partition.size(); i++) {
    current_partition[i].push_back(set[index]);
    get_partitions(set, index + 1, ans, current_partition);
    current_partition[i].pop_back();
  }

  // Add the current element as a
  // singleton subset and recall
  current_partition.push_back({set[index]});
  get_partitions(set, index + 1, ans, current_partition);
  current_partition.pop_back();
}

// Helper function to print the results
void print_partitions(const all_partitions &partitions) {
  for (const auto &p : partitions) {
    std::cout << "{ ";
    for (const auto &s : p) {
      std::cout << "{ ";
      for (int val : s) { std::cout << val << " "; }
      std::cout << "} ";
    }
    std::cout << "}\n";
  }
}

std::vector<std::vector<int>> generate_permutations(std::vector<int> &set) {

  std::vector<std::vector<int>> permutations;
  std::sort(set.begin(), set.end());
  do { permutations.push_back(set); } while (std::next_permutation(set.begin(), set.end()));
  return permutations;
}

int permutation_sign(std::vector<int> perm) {

  int inversions = 0;
  for (size_t i = 0; i < perm.size(); ++i) {
    for (size_t j = i + 1; j < perm.size(); ++j) {
      if (perm[i] > perm[j]) { ++inversions; }
    }
  }
  return (inversions % 2 == 0) ? 1 : -1;
}

void cumulant_decomposition(std::vector<int> &unprimed, std::vector<int> &primed) {
  all_partitions partitions;
  partition c_partition;
  get_partitions(unprimed, 0, partitions, c_partition);

  auto primed_perms = generate_permutations(primed);

  std::cout << "G(";
  for (int val : unprimed) { std::cout << val; }
  std::cout << "|";
  for (int val : primed) { std::cout << val << "'"; }
  std::cout << ")" << std::endl;

  for (auto current_partition : partitions) {
    if (current_partition.size() == 1) continue; // Skip the full partition

    std::set<std::vector<int>> visited_permutations;
    for (auto perm : primed_perms) {

      int primed_cursor = 0;
      std::vector<int> effective_permutation;
      std::vector<subset> term_subsets; // Store the generated subsets for printing

      // 1. Build the canonical permutation and the term's subsets in one go.
      for (auto &x : current_partition) {
        std::vector<int> primed_perm_slice(perm.begin() + primed_cursor, perm.begin() + primed_cursor + x.size());
        primed_cursor += x.size();

        std::sort(primed_perm_slice.begin(), primed_perm_slice.end()); // Canonical form

        // Store the sorted slice for printing later
        term_subsets.push_back(primed_perm_slice);

        // Add the elements to the effective_permutation for the uniqueness check
        effective_permutation.insert(effective_permutation.end(), primed_perm_slice.begin(), primed_perm_slice.end());
      }

      // 2. Check if this canonical permutation has already been processed.
      if (visited_permutations.count(effective_permutation)) { continue; }
      visited_permutations.insert(effective_permutation);

      // 3. If it's new, calculate its sign and print the full term.
      int sign = permutation_sign(effective_permutation);
      std::cout << (sign == 1 ? "- " : "+"); //opposite sign convention

      for (size_t i = 0; i < current_partition.size(); ++i) {
        const auto &unprimed_subset = current_partition[i];
        const auto &primed_subset   = term_subsets[i]; // Use the stored subset

        std::cout << "C(";
        for (int val : unprimed_subset) { std::cout << val; }
        std::cout << "|";
        for (int val : primed_subset) { std::cout << val << "'"; }
        std::cout << ") ";
      }
      std::cout << std::endl;
    }
  }
}

int main() {
  std::vector<int> unprimed = {1, 2, 3};
  std::vector<int> primed   = {1, 2, 3};

  cumulant_decomposition(unprimed, primed);
  return 0;
}