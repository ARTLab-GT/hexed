#include "catch.hpp"

#include "Grid.hpp"
#include "Equidistant.hpp"

TEST_CASE("Grid")
{
  Equidistant basis (8);
  Grid grid1 (4, 1, 5, 0.1, basis);
  Grid grid2 (4, 2, 5, 0.1, basis);
  Grid grid3 (4, 3, 27, 0.1, basis);
  SECTION("Construction")
  {
    REQUIRE(grid2.basis.node(1) == 1./7.);
    REQUIRE(grid3.n_var == 4);
    REQUIRE(grid1.n_dim == 1);
    REQUIRE(grid2.n_dim == 2);
    REQUIRE(grid3.n_dim == 3);
    REQUIRE(grid1.n_qpoint == 8);
    REQUIRE(grid2.n_qpoint == 64);
    REQUIRE(grid3.n_qpoint == 512);
    REQUIRE(grid1.n_dof == 32);
    REQUIRE(grid2.n_dof == 256);
    REQUIRE(grid3.n_dof == 2048);
    REQUIRE(grid1.n_elem == 5);
    int size = 2048*27;
    REQUIRE(grid3.state_r[0] == 0.);
    REQUIRE(grid3.state_w[0] == 0.);
    REQUIRE(grid3.state_r[size - 1] == 0.);
    REQUIRE(grid3.state_w[size - 1] == 0.);
  }
  SECTION("Positioning and visualization")
  {
  }
}
