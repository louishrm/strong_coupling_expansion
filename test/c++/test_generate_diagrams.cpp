#include <gtest/gtest.h>
#include "sc_expansion/generate_diagrams.hpp"
#include <vector>
#include <iostream>

using namespace sc_expansion;

TEST(GenerateDiagramsTest, Order2Diagrams) {
  VacuumDiagramGenerator gen(2);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  int count_nonzero_fm = 0;
  for (const auto &g : unique_graphs) {
    int V = std::sqrt(g.size());
    Graph graph(g, V);
    if (graph.get_free_multiplicity() > 0) count_nonzero_fm++;
  }

  EXPECT_EQ(count_nonzero_fm, 1);
}

TEST(GenerateDiagramsTest, Order4Diagrams) {
  VacuumDiagramGenerator gen(4);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  int count_nonzero_fm = 0;
  for (const auto &g : unique_graphs) {
    int V = std::sqrt(g.size());
    Graph graph(g, V);
    if (graph.get_free_multiplicity() > 0) count_nonzero_fm++;
  }

  EXPECT_EQ(count_nonzero_fm, 3);
}

TEST(GenerateDiagramsTest, Order6Diagrams) {
  VacuumDiagramGenerator gen(6);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  int count_nonzero_fm = 0;
  for (const auto &g : unique_graphs) {
    int V = std::sqrt(g.size());
    Graph graph(g, V);
    if (graph.get_free_multiplicity() > 0) count_nonzero_fm++;
  }

  EXPECT_EQ(count_nonzero_fm, 7);
  std::cout << "Detected " << count_nonzero_fm << " unique diagrams of order 6 with non-zero FM." << std::endl;
}

TEST(GenerateDiagramsTest, Order8Diagrams) {
  VacuumDiagramGenerator gen(8);
  gen.generate();

  const auto &unique_graphs = gen.get_unique_graphs();
  int count_nonzero_fm = 0;
  for (const auto &g : unique_graphs) {
    int V = std::sqrt(g.size());
    Graph graph(g, V);
    if (graph.get_free_multiplicity() > 0) count_nonzero_fm++;
  }

  std::cout << "Detected " << count_nonzero_fm << " unique diagrams of order 8 with non-zero FM out of " << unique_graphs.size() << " unique topologies." << std::endl;
}
