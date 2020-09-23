#include <catch.hpp>

#include "Solution.hpp"

TEST_CASE("Solution class")
{
  Solution sol (3, 2, 4, 0.7);
  Grid* g;
  REQUIRE_THROWS(g = &sol.get_grid(0));
  std::vector<int> lc {-1, -1};
  std::vector<int> uc {1, 2};
  sol.add_block_grid(1, lc, uc);
  g = &sol.get_grid(0);
  REQUIRE(g->basis.rank == 4);
  REQUIRE(g->n_elem == 6);
}
