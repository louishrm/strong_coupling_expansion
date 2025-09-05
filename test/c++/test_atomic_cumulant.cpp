#include "../c++/sc_expansion/atomic_cumulant.hpp"
#include <vector>
#include <iostream>
#include <string>

int main(int argc, char *argv[]) {

  std::vector<std::string> unprimed = {"1", "2", "3"};
  std::vector<std::string> primed   = {"1'", "2'", "3'"};

  atomic_cumulant::showCumulantDecomposition(unprimed, primed);

  // // Generate all partitions of the 'unprimed' vector.
  // atomic_cumulant::AllPartitions partitions = atomic_cumulant::allPartitions(unprimed);

  // // Print a header for the output.
  // std::cout << "All partitions of {" << atomic_cumulant::format_vector(unprimed) << "}:\n";

  // // Iterate through each partition in the result.
  // for (const auto &partition : partitions) {
  //   std::cout << "{ ";
  //   // Iterate through each subset within the partition.
  //   for (size_t i = 0; i < partition.size(); ++i) {
  //     const auto &subset = partition[i];
  //     // Use format_vector to print the contents of the subset.
  //     std::cout << "{" << atomic_cumulant::format_vector(subset) << "}";
  //     // Add a comma between subsets for readability.
  //     if (i < partition.size() - 1) { std::cout << ", "; }
  //   }
  //   std::cout << " }\n";

  return 0;
}