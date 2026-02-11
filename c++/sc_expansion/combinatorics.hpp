#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdint>

namespace sc_expansion {

  inline uint64_t factorial(uint64_t k) {
    if (k <= 1) return 1;
    return k * factorial(k - 1);
  }

  class PartitionGenerator {
    public:
    //calculates all valid partitions immediately
    PartitionGenerator(int n, int max_val) : current_index(0) {
      std::vector<int> buffer(n); // Buffer to hold work-in-progress
      generate_recursive(n, buffer, 0, max_val);
    }

    // Access the current partition
    const std::vector<int> &current() const { return all_partitions[current_index]; }

    // Advance to the next partition. Returns false if we are at the end.
    bool next() {
      current_index++;
      return current_index < all_partitions.size();
    }

    // Check if valid (useful for initial check)
    bool is_valid() const { return current_index < all_partitions.size(); }

    // Reset iterator to reuse the generator
    void reset() { current_index = 0; }

    private:
    std::vector<std::vector<int>> all_partitions;
    size_t current_index;

    // Your original recursive logic, adapted to store results
    void generate_recursive(int n, std::vector<int> &v, int level, int max_val) {
      if (n == 0) {
        // BASE CASE: We found a valid partition.
        // Save exactly the used portion of the buffer (0 to level)
        if (level > 0) {
          std::vector<int> valid_part(v.begin(), v.begin() + level);
          all_partitions.push_back(valid_part);
        }
        return;
      }

      // Determine starting number to ensure non-decreasing order
      int first = (level == 0) ? 1 : v[level - 1];

      for (int i = first; i <= n && i <= max_val; i++) {
        v[level] = i;
        generate_recursive(n - i, v, level + 1, max_val);
      }
    }
  };

} // namespace sc_expansion