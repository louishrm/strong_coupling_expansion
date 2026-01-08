#include "cumulant.hpp"
#include <numeric>
#include <set>
#include <algorithm>

namespace {

  // Type definitions for clarity
  using subset_t         = std::vector<int>;
  using partition_t      = std::vector<subset_t>;
  using all_partitions_t = std::vector<partition_t>;

  void fill_partitions(std::vector<int> const &set, int index, all_partitions_t &ans, partition_t &current_partition) {

    if (index == set.size()) {
      ans.push_back(current_partition);
      return;
    }

    // Option 1: Add the current element to an existing subset
    for (auto &subset : std::vector<subset_t>(current_partition)) {
      subset.push_back(set[index]);
      fill_partitions(set, index + 1, ans, current_partition);
      subset.pop_back(); // Backtrack
    }

    // Option 2: Add the current element as a new singleton subset
    current_partition.push_back({set[index]});
    fill_partitions(set, index + 1, ans, current_partition);
    current_partition.pop_back(); // Backtrack
  }

  all_partitions_t get_partitions(std::vector<int> const &set) {

    all_partitions_t ans;
    partition_t current_partition;
    fill_partitions(set, 0, ans, current_partition);
    return ans;
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

  bool is_spin_conserving(hubbard_atom::cumul_args const &unprimed, hubbard_atom::cumul_args const &primed) {
    int up_count_unprimed   = 0;
    int down_count_unprimed = 0;
    for (const auto &arg : unprimed) {
      if (arg.second == 0)
        up_count_unprimed++;
      else
        down_count_unprimed++;
    }

    int up_count_primed   = 0;
    int down_count_primed = 0;
    for (const auto &arg : primed) {
      if (arg.second == 0)
        up_count_primed++;
      else
        down_count_primed++;
    }

    return (up_count_unprimed == up_count_primed) && (down_count_unprimed == down_count_primed);
  }

} // namespace

namespace sc_expansion {

  double compute_cumulant_decomposition(hubbard_atom::cumul_args const &unprimed, hubbard_atom::cumul_args const &primed,
                                        triqs::atom_diag::atom_diag<false> const &ad, double beta) {

    /*Core logic: 
  -This is a recursive function that should compute C^0_n( unprimed | primed).
  1. Create index set = {0, 1, 2, ..., n-1} where n = size of unprimed (or primed, they should be equal).
  2. First we check if unprimed and primed have size 1. Then we return G^0_1( unprimed | primed).
  3. The first term is G^0_n( unprimed |primed) which we can compute directly.
  4. Initialize a result variable that computes C^0_n - G^0_n.
  5. Generate all partitions of the unprimed indices, and all the permutations of the primed indices. 
  6. Loop through partitions and permutations (skip the trivial partition). 
  7. Make a set to track visited canonical permutations.
  8. For each permutation, create an array to store the resulting "effective permutation" and a cursor to track the subset of the unprimed partition and an array to store the primed subsets. 
  9. Loop through each subset of the current partition. For each subset, take a slice of the primed permutation corresponding to the size of the current unprimed subset. Sort this slice to get the canonical form. Add these elements to the effective permutation array and store the sorted slice in the primed subsets array.
  10. After constructing the effective permutation, check if it has already been visited. If so, skip to the next permutation. If not, mark it as visited.
  11. Calculate the sign and initialize the corresponding products contrbution. 
  12. Build the argument lists for the recursive calls and recursively call the function for each subset pair. Multiply the results to get the current product.
  13. Add the signed product to the result variable.
  14. Finally return G^0_n + result.

  TODO: Optimize this function to cache results for specific (unprimed, primed) pairs to avoid redundant calculations.

  */

    bool spin_cons = is_spin_conserving(unprimed, primed);
    if (!spin_cons) { return 0.0; }

    //1 Create index vectors for unprimed and primed args.
    int order = unprimed.size();
    std::vector<int> unprimed_indices(order);
    std::iota(unprimed_indices.begin(), unprimed_indices.end(), 0);

    std::vector<int> primed_indices(order);
    std::iota(primed_indices.begin(), primed_indices.end(), 0);

    //2 Base case: if size is 1, return G0
    if (unprimed.size() == 1) { return hubbard_atom::G0(ad, beta, unprimed, primed); }

    //3 Compute the first term: G^0_n( unprimed | primed)
    double G0n = hubbard_atom::G0(ad, beta, unprimed, primed);

    //4 initializae the subtraction term.
    double low_order_cumulants = 0.0;

    //5 Generate all partitions of unprimed indices and permutations of primed indices.
    all_partitions_t unprimed_partitions = get_partitions(unprimed_indices);

    auto primed_perms = generate_permutations(primed_indices);

    //6 Loop through partitions and permutations.
    for (partition_t current_partition : unprimed_partitions) {
      if (current_partition.size() == 1) continue; // Skip the trivial partition.

      //7
      std::set<std::vector<int>> visited_permutations;

      for (auto perm : primed_perms) {

        //8
        int primed_cursor = 0;
        std::vector<int> effective_permutation(order);
        std::vector<subset_t> primed_subsets;
        std::vector<int>::iterator write_cursor = effective_permutation.begin();

        //9.
        for (auto x : current_partition) {

          auto write_end_cursor = write_cursor + x.size();
          std::copy(perm.begin() + primed_cursor, perm.begin() + primed_cursor + x.size(), write_cursor);
          std::sort(write_cursor, write_end_cursor);
          primed_subsets.push_back(std::vector<int>(write_cursor, write_end_cursor));
          primed_cursor += x.size();
          write_cursor = write_end_cursor;
        }

        //10
        if (visited_permutations.count(effective_permutation)) { continue; }
        visited_permutations.insert(effective_permutation);

        //11
        double current_product = 1.0;
        int sign               = -permutation_sign(effective_permutation); //opposite sign convention

        //12
        for (size_t i = 0; i < current_partition.size(); ++i) {
          const auto &unprimed_idx_subset = current_partition[i];
          const auto &primed_idx_subset   = primed_subsets[i];

          hubbard_atom::cumul_args current_unprimed_args;
          current_unprimed_args.reserve(unprimed_idx_subset.size());
          for (int index : unprimed_idx_subset) { current_unprimed_args.push_back(unprimed[index]); }

          hubbard_atom::cumul_args current_primed_args;
          current_primed_args.reserve(primed_idx_subset.size());
          for (int index : primed_idx_subset) { current_primed_args.push_back(primed[index]); }

          current_product *= compute_cumulant_decomposition(current_unprimed_args, current_primed_args, ad, beta);
        }
        low_order_cumulants += sign * current_product;
      }
    }
    return G0n + low_order_cumulants;
  }

  //COMMENTED OUT FOR NOW, add to tests later

  // void cumulant_decomposition(std::vector<int> &unprimed, std::vector<int> &primed) {
  //   all_partitions partitions;
  //   partition c_partition;
  //   get_partitions(unprimed, 0, partitions, c_partition);

  //   auto primed_perms = generate_permutations(primed);

  //   std::cout << "G(";
  //   for (int val : unprimed) { std::cout << val; }
  //   std::cout << "|";
  //   for (int val : primed) { std::cout << val << "'"; }
  //   std::cout << ")" << std::endl;

  //   for (auto current_partition : partitions) {
  //     if (current_partition.size() == 1) continue; // Skip the full partition

  //     std::set<std::vector<int>> visited_permutations;
  //     for (auto perm : primed_perms) {

  //       int primed_cursor = 0;
  //       std::vector<int> effective_permutation;
  //       std::vector<subset> term_subsets; // Store the generated subsets for printing

  //       // 1. Build the canonical permutation and the term's subsets in one go.
  //       for (auto &x : current_partition) {
  //         std::vector<int> primed_perm_slice(perm.begin() + primed_cursor, perm.begin() + primed_cursor + x.size());
  //         primed_cursor += x.size();

  //         std::sort(primed_perm_slice.begin(), primed_perm_slice.end()); // Canonical form

  //         // Store the sorted slice for printing later
  //         term_subsets.push_back(primed_perm_slice);

  //         // Add the elements to the effective_permutation for the uniqueness check
  //         effective_permutation.insert(effective_permutation.end(), primed_perm_slice.begin(), primed_perm_slice.end());
  //       }

  //       // 2. Check if this canonical permutation has already been processed.
  //       if (visited_permutations.count(effective_permutation)) { continue; }
  //       visited_permutations.insert(effective_permutation);

  //       // 3. If it's new, calculate its sign and print the full term.
  //       int sign = permutation_sign(effective_permutation);
  //       std::cout << (sign == 1 ? "- " : "+"); //opposite sign convention

  //       for (size_t i = 0; i < current_partition.size(); ++i) {
  //         const auto &unprimed_subset = current_partition[i];
  //         const auto &primed_subset   = term_subsets[i]; // Use the stored subset

  //         std::cout << "C(";
  //         for (int val : unprimed_subset) { std::cout << val; }
  //         std::cout << "|";
  //         for (int val : primed_subset) { std::cout << val << "'"; }
  //         std::cout << ") ";
  //       }
  //       std::cout << std::endl;
  //     }
  //   }
  // }

  // // Helper function to print the results
  // void print_partitions(const all_partitions &partitions) {
  //   for (const auto &p : partitions) {
  //     std::cout << "{ ";
  //     for (const auto &s : p) {
  //       std::cout << "{ ";
  //       for (int val : s) { std::cout << val << " "; }
  //       std::cout << "} ";
  //     }
  //     std::cout << "}\n";
  //   }
  // }

} // namespace sc_expansion
