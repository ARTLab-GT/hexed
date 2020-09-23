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
  sol.add_block_grid(2);
  g = &sol.get_grid(0);
  REQUIRE(g->basis.rank == 4);
  REQUIRE(g->n_elem == 6);
  REQUIRE(g->mesh_size == 0.35);
  g = &sol.get_grid(1);
  REQUIRE(g->basis.rank == 4);
  REQUIRE(g->n_elem == 16);
  REQUIRE(g->mesh_size == 0.175);
  REQUIRE(g->pos[0] == 0);
}
