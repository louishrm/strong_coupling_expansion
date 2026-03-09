#pragma once
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdint>

namespace sc_expansion {

  inline double compute_permutation_sign(std::vector<int> const &p) {
    int n            = p.size();
    int cycles       = 0;
    uint64_t visited = 0; // Each bit represents an index (up to 64 operators)

    for (int i = 0; i < n; ++i) {
      if (!(visited & (1ULL << i))) { // Check if the i-th bit is 0
        int current = i;
        while (!(visited & (1ULL << current))) {
          visited |= (1ULL << current); // Set the current bit to 1
          current = p[current];
        }
        cycles++;
      }
    }
    return ((n - cycles) & 1) ? -1 : 1;
  }

  inline uint64_t factorial(uint64_t k) {
    if (k <= 1) return 1;
    return k * factorial(k - 1);
  }

  class SJT {

    public:
    SJT(int n_) : n(n_) {

      this->permutation.resize(n);
      std::iota(this->permutation.begin(), this->permutation.end(), 1);

      this->directions.resize(n + 1, -1); // -1 for left, +1 for right //use index from 0 to n included

      this->positions.resize(n + 1);
      for (int v = 1; v <= n; v++) { this->positions[v] = v - 1; }
    }

    const std::vector<int> &get_permutation() const { return this->permutation; }

    bool next_permutation() {

      /*Find the largest element with non zero direction
      Move it in the direction specified by the directions vector
      */
      for (int i = n; i >= 1; i--) {
        if (this->is_mobile(i)) {
          this->perform_swap(i);
          this->reverse_directions(i);
          return true;
        }
      }
      return false;
    }

    private:
    int n;
    std::vector<int> permutation;
    std::vector<int> directions;
    std::vector<int> positions;

    bool is_mobile(int i) const {

      /*is mobile if direction is not zero 
      and can move  within the bounds without hitting a larger number*/
      int dir        = this->directions[i];
      int idx        = this->positions[i];
      int target_idx = idx + dir;
      if (idx + dir < 0 || idx + dir >= n) { return false; }
      return this->permutation[target_idx] < this->permutation[idx];
    }

    void perform_swap(int i) {

      int idx        = this->positions[i];
      int target_idx = idx + this->directions[i];
      // Swap in the permutation
      std::swap(this->permutation[idx], this->permutation[target_idx]);
      //update positions
      this->positions[i]                      = target_idx;
      this->positions[this->permutation[idx]] = idx;
    }

    void reverse_directions(int i) {
      for (int j = i + 1; j <= n; j++) { this->directions[j] *= -1; }
    }
  };

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