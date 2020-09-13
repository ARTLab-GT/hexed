#include "catch.hpp"

#include "Grid.hpp"
#include "Equidistant.hpp"

TEST_CASE("Grid")
{
  Equidistant basis (8);
  Grid grid (4, 3, 5, basis);
  REQUIRE(grid.basis.node(1) == 1./7.);
  REQUIRE(grid.n_var == 4);
  REQUIRE(grid.n_dim == 3);
  REQUIRE(grid.n_qpoint == 512);
  REQUIRE(grid.n_dof == 2048);
  REQUIRE(grid.n_elem == 5);
  int size = 2048*5;
  REQUIRE(grid.access_r[0] == 0.);
  REQUIRE(grid.access_w[0] == 0.);
  REQUIRE(grid.access_r[size - 1] == 0.);
  REQUIRE(grid.access_w[size - 1] == 0.);
  grid.access_r[7] = 2.5;
  REQUIRE(grid.access_r[7] == 2.5);
}
